import numpy as np

from cndpy.transform import compute_vx
from cndpy.variance import cumulative_series, get_cumulative_variance


def test_iteration_count_is_n_minus_3(nopal_df, nopal_nutrient_cols):
    df, _ = compute_vx(nopal_df, nopal_nutrient_cols)
    n = len(df)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    n_vcols = len(nopal_nutrient_cols) + 1  # + V_R
    expected_iterations = n - 3  # range(2, n-1) has n-3 values
    assert len(var_df) == expected_iterations * n_vcols


def test_ratio_is_never_negative(nopal_df, nopal_nutrient_cols):
    df, _ = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    assert (var_df['f_i'] >= 0).all()


def test_zero_high_variance_gives_zero_ratio():
    # Construct a case where the HIGH group's V_A is identical for the first
    # two (highest-yield) rows: at i=2, var_n1 == 0 (ddof=1 var of 2 equal
    # values), which the code guards against and forces ratio == 0.0.
    import pandas as pd

    df = pd.DataFrame({
        'yield': [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0],
        'A': [5.0, 5.0, 4.0, 3.0, 2.0, 6.0, 7.0, 8.0],
    })
    df, _ = compute_vx(df, ['A'])
    var_df = get_cumulative_variance(df, ['A'])
    first_row = var_df.iloc[0]
    assert first_row['f_i'] == 0.0


def test_cumulative_series_ends_at_100(nopal_df, nopal_nutrient_cols):
    df, _ = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    sub = cumulative_series(var_df, 'V_N')
    assert np.isclose(sub['FiC'].iloc[-1], 100.0)


def test_cumulative_series_flat_zero_when_total_nonpositive():
    import pandas as pd

    var_df = pd.DataFrame({
        'yield_cut': [10, 9, 8],
        'nut': ['V_X', 'V_X', 'V_X'],
        'f_i': [0.0, 0.0, 0.0],
    })
    sub = cumulative_series(var_df, 'V_X')
    assert (sub['FiC'] == 0.0).all()


def test_cumulative_series_is_monotonically_nondecreasing(nopal_df, nopal_nutrient_cols):
    df, _ = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    sub = cumulative_series(var_df, 'V_N')
    diffs = np.diff(sub['FiC'].values)
    assert (diffs >= -1e-9).all()
