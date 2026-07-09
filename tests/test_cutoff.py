import numpy as np
import pandas as pd
import pytest

from cndpy.cutoff import fit_cubic_inflections, select_yield_cutoff
from cndpy.transform import compute_vx
from cndpy.variance import get_cumulative_variance


def test_fit_cubic_inflections_skips_columns_with_few_points():
    var_df = pd.DataFrame({
        'yield_cut': [10, 9, 8],
        'nut': ['V_A', 'V_A', 'V_A'],
        'f_i': [1.0, 1.0, 1.0],
    })
    fits = fit_cubic_inflections(var_df, ['A'])
    assert 'A' not in fits


def test_select_yield_cutoff_recovers_known_inflection_point():
    # FiC(Y) = Y^3 - 15Y^2 -> Y* = -b/(3a) = 15/3 = 5
    coeffs = np.array([1.0, -15.0, 0.0, 0.0])
    x_vals = np.array([0.0, 10.0])
    y_vals = np.polyval(coeffs, x_vals)
    cubic_fits = {'X': (coeffs, 0.99, x_vals, y_vals)}

    df = pd.DataFrame({'yield': [10, 8, 6, 4, 2, 1]})
    result = select_yield_cutoff(df, cubic_fits)

    assert result['y_star_ref'] == pytest.approx(5.0)
    assert result['used_fallback'] is False
    assert result['best_n1'] == 4
    assert result['cutoff'] == pytest.approx(4.0)
    assert list(result['high_df']['yield']) == [10, 8, 6, 4]


def test_select_yield_cutoff_falls_back_to_75th_percentile_with_no_candidates():
    df = pd.DataFrame({'yield': [10, 8, 6, 4, 2, 1]})
    result = select_yield_cutoff(df, cubic_fits={})
    assert result['used_fallback'] is True
    assert result['warning'] is not None
    assert result['y_star_ref'] == pytest.approx(df['yield'].quantile(0.75))


def test_nopal_regression_cutoff_and_n1(nopal_df, nopal_nutrient_cols):
    """Regression test against the 2026-07-08 verification of CND_V14.py
    against TODO-nopal.xls: cutoff=34.2 kg/plant, n1=14 (38.9%)."""
    df, _d_plus_1 = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    cubic_fits = fit_cubic_inflections(var_df, nopal_nutrient_cols)
    result = select_yield_cutoff(df, cubic_fits)

    assert result['cutoff'] == pytest.approx(34.2, abs=0.05)
    assert result['best_n1'] == 14

    # N, Ca, Mg, R are the documented in-context candidates; P and K are
    # out of context (Y* far outside the observed yield range) and must be
    # absent from cubic_fits' contribution to the mean, but they may still
    # appear in cubic_fits itself (fit_cubic_inflections doesn't filter by
    # range — only select_yield_cutoff does).
    for col in ('N', 'Ca', 'Mg', 'R'):
        assert col in cubic_fits
