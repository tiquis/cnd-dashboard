import pytest

from cndpy.cutoff import fit_cubic_inflections, select_yield_cutoff
from cndpy.norms import compute_norms, critical_chi_square, low_yield_proportion
from cndpy.transform import compute_vx
from cndpy.variance import get_cumulative_variance

# Magallanes-Quintanar et al. (2004), published V*_X norms for nopal (n=36).
PUBLISHED_NOPAL_NORMS = {
    'V_N': -1.1487, 'V_P': -2.25402, 'V_K': 0.36933,
    'V_Ca': 0.3696, 'V_Mg': -0.72992, 'V_R': 3.39372,
}


def test_nopal_norms_match_published_values(nopal_df, nopal_nutrient_cols):
    df, _d_plus_1 = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    cubic_fits = fit_cubic_inflections(var_df, nopal_nutrient_cols)
    result = select_yield_cutoff(df, cubic_fits)
    norms = compute_norms(result['high_df'], nopal_nutrient_cols)

    for key, published in PUBLISHED_NOPAL_NORMS.items():
        assert norms[key] == pytest.approx(published, abs=0.0153)


def test_nopal_norms_sum_to_zero(nopal_df, nopal_nutrient_cols):
    df, _d_plus_1 = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    cubic_fits = fit_cubic_inflections(var_df, nopal_nutrient_cols)
    result = select_yield_cutoff(df, cubic_fits)
    norms = compute_norms(result['high_df'], nopal_nutrient_cols)

    v_sum = sum(norms[f'V_{c}'] for c in nopal_nutrient_cols + ['R'])
    assert v_sum == pytest.approx(0.0, abs=1e-3)


def test_critical_chi_square_matches_nopal_reference():
    assert critical_chi_square(prop_low=0.611, d_plus_1=6) == pytest.approx(4.4867, abs=0.001)


def test_low_yield_proportion():
    assert low_yield_proportion(n_total=36, n_high=14) == pytest.approx(1 - 14 / 36)
