"""End-to-end safety net: reproduces the values from the 2026-07-08
verification of CND_V14.py against the original TODO-nopal.xls worksheet
(see cndpy_refactor_plan.md section 0 / Resumen_Proyecto_CND_Dashboard.md).
If a future change to cndpy breaks the pipeline, this file must fail first.
"""
import pytest

from cndpy.cutoff import fit_cubic_inflections, select_yield_cutoff
from cndpy.diagnosis import compute_indices
from cndpy.norms import (
    chisquare_goodness_of_fit,
    compute_norms,
    critical_chi_square,
    low_yield_proportion,
)
from cndpy.transform import compute_vx
from cndpy.variance import get_cumulative_variance

PUBLISHED_NOPAL_NORMS = {
    'V_N': -1.1487, 'V_P': -2.25402, 'V_K': 0.36933,
    'V_Ca': 0.3696, 'V_Mg': -0.72992, 'V_R': 3.39372,
}


def _run_pipeline(df_raw, nutrient_cols):
    df, d_plus_1 = compute_vx(df_raw, nutrient_cols)
    var_df = get_cumulative_variance(df, nutrient_cols)
    cubic_fits = fit_cubic_inflections(var_df, nutrient_cols)
    cutoff_result = select_yield_cutoff(df, cubic_fits)
    norms = compute_norms(cutoff_result['high_df'], nutrient_cols)
    prop_low = low_yield_proportion(len(df), len(cutoff_result['high_df']))
    critical_r2 = critical_chi_square(prop_low, d_plus_1)
    indices_df = compute_indices(df, nutrient_cols, norms, critical_r2)
    return dict(
        df=df, d_plus_1=d_plus_1, cubic_fits=cubic_fits, cutoff_result=cutoff_result,
        norms=norms, prop_low=prop_low, critical_r2=critical_r2, indices_df=indices_df,
    )


def test_nopal_end_to_end_reproduces_2026_07_08_verification(nopal_df, nopal_nutrient_cols):
    ctx = _run_pipeline(nopal_df, nopal_nutrient_cols)

    assert ctx['cutoff_result']['best_n1'] == 14
    assert ctx['cutoff_result']['cutoff'] == pytest.approx(34.2, abs=0.05)
    assert ctx['prop_low'] == pytest.approx(0.611, abs=0.001)
    assert ctx['critical_r2'] == pytest.approx(4.4867, abs=0.001)
    assert ctx['d_plus_1'] == 6

    for key, published in PUBLISHED_NOPAL_NORMS.items():
        assert ctx['norms'][key] == pytest.approx(published, abs=0.0153)

    r2_fit = chisquare_goodness_of_fit(ctx['indices_df'], ctx['d_plus_1'])
    assert r2_fit == pytest.approx(0.9520, abs=0.001)

    # Documented in-context candidates for the cutoff mean: N, Ca, Mg, R
    # (P and K are out of context, Y*≈250.6 and 78.4 kg/plant respectively).
    for col in ('N', 'Ca', 'Mg', 'R'):
        assert col in ctx['cubic_fits']


def test_maiz_end_to_end_pipeline_is_self_consistent(maiz_df, maiz_nutrient_cols):
    """Sanity/self-consistency check for the maize dataset.

    NOTE: unlike nopal, no literal published V*_X table for maize was
    available to re-verify in this session (the plan documents differences
    of 0.0010-0.0140 vs. Magallanes-Quintanar et al. 2006, but that
    cross-check was not independently redone here - see plan section 0).
    This test therefore only asserts internal consistency, not agreement
    with the published manuscript values.
    """
    ctx = _run_pipeline(maiz_df, maiz_nutrient_cols)

    assert ctx['d_plus_1'] == len(maiz_nutrient_cols) + 1
    assert 0 < ctx['cutoff_result']['best_n1'] < len(maiz_df)
    v_sum = sum(ctx['norms'][f'V_{c}'] for c in maiz_nutrient_cols + ['R'])
    assert v_sum == pytest.approx(0.0, abs=1e-3)
    assert ctx['cutoff_result']['cutoff'] == pytest.approx(
        ctx['cutoff_result']['high_df']['yield'].min()
    )
