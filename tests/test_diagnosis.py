import pandas as pd
import pytest

from cndpy.diagnosis import compute_indices, most_limiting_nutrient, rank_limiting_nutrients

NUTRIENT_COLS = ['A', 'B']
NORMS = {
    'V_A': 0.0, 'SD_A': 1.0,
    'V_B': 0.5, 'SD_B': 2.0,
    'V_R': -0.5, 'SD_R': 0.5,
}


def _make_vx_df():
    return pd.DataFrame({
        'V_A': [1.0, -2.0],
        'V_B': [0.5, 3.0],
        'V_R': [-1.0, 0.0],
    })


def test_compute_indices_row_matches_batch_vs_single():
    df_vx = _make_vx_df()
    batch = compute_indices(df_vx, NUTRIENT_COLS, NORMS, critical_r2=5.0)

    for i in range(len(df_vx)):
        single = compute_indices(df_vx.iloc[[i]], NUTRIENT_COLS, NORMS, critical_r2=5.0)
        for col in ['I_A', 'I_B', 'I_R', 'CND_r2', 'Balanced']:
            assert single[col].iloc[0] == pytest.approx(batch[col].iloc[i]) \
                if col != 'Balanced' else single[col].iloc[0] == batch[col].iloc[i]


def test_compute_indices_values_and_r2():
    df_vx = _make_vx_df()
    result = compute_indices(df_vx, NUTRIENT_COLS, NORMS, critical_r2=100.0)

    # Row 0: I_A=(1.0-0.0)/1.0=1.0, I_B=(0.5-0.5)/2.0=0.0, I_R=(-1.0-(-0.5))/0.5=-1.0
    assert result['I_A'].iloc[0] == pytest.approx(1.0)
    assert result['I_B'].iloc[0] == pytest.approx(0.0)
    assert result['I_R'].iloc[0] == pytest.approx(-1.0)
    assert result['CND_r2'].iloc[0] == pytest.approx(1.0 ** 2 + 0.0 ** 2 + 1.0 ** 2)
    assert result['Balanced'].iloc[0] == (result['CND_r2'].iloc[0] < 100.0)


def test_rank_limiting_nutrients_orders_by_descending_abs_value():
    row = pd.Series({'I_A': 0.1, 'I_B': -5.0, 'I_R': 2.0})
    ranked = rank_limiting_nutrients(row, NUTRIENT_COLS)
    assert ranked == 'B > Rd > A'


def test_most_limiting_nutrient_identifies_max_abs_column():
    df = pd.DataFrame({'I_A': [0.1, 9.0], 'I_B': [-5.0, 0.2], 'I_R': [2.0, -0.5]})
    result = most_limiting_nutrient(df, NUTRIENT_COLS)
    assert list(result) == ['B', 'A']
