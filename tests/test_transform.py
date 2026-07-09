import math

import numpy as np
import pandas as pd

from cndpy.transform import compute_vx


def test_sum_of_vx_is_zero_per_row(nopal_df, nopal_nutrient_cols):
    df, d_plus_1 = compute_vx(nopal_df, nopal_nutrient_cols)
    v_cols = [f'V_{c}' for c in nopal_nutrient_cols] + ['V_R']
    row_sums = df[v_cols].sum(axis=1)
    assert np.allclose(row_sums, 0.0, atol=1e-9)


def test_filling_value_is_complement_to_100(nopal_df, nopal_nutrient_cols):
    df, _ = compute_vx(nopal_df, nopal_nutrient_cols)
    assert np.allclose(df['R'], 100 - df[nopal_nutrient_cols].sum(axis=1))


def test_d_plus_1_equals_nutrient_count_plus_one(nopal_df, nopal_nutrient_cols):
    _, d_plus_1 = compute_vx(nopal_df, nopal_nutrient_cols)
    assert d_plus_1 == len(nopal_nutrient_cols) + 1


def test_known_two_nutrient_case():
    # A=20, B=30 -> R = 100-50 = 50; components [20, 30, 50]
    # g = (20*30*50)**(1/3), computed independently via math.log/exp
    df = pd.DataFrame({'A': [20.0], 'B': [30.0]})
    out, d_plus_1 = compute_vx(df, ['A', 'B'])

    assert d_plus_1 == 3
    g = (20.0 * 30.0 * 50.0) ** (1.0 / 3.0)
    expected_v_a = math.log(20.0 / g)
    expected_v_b = math.log(30.0 / g)
    expected_v_r = math.log(50.0 / g)

    assert out['R'].iloc[0] == 50.0
    assert math.isclose(out['V_A'].iloc[0], expected_v_a, rel_tol=1e-12)
    assert math.isclose(out['V_B'].iloc[0], expected_v_b, rel_tol=1e-12)
    assert math.isclose(out['V_R'].iloc[0], expected_v_r, rel_tol=1e-12)
