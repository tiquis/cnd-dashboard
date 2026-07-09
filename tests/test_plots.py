import pandas as pd
import plotly.graph_objects as go

from cndpy.cutoff import (
    fit_cubic_inflections,
    select_representative_column,
    select_yield_cutoff,
)
from cndpy.diagnosis import compute_indices
from cndpy.norms import compute_norms, critical_chi_square, low_yield_proportion
from cndpy.plots import chisquare_cdf_figure, cumulative_variance_figure, diagnosis_bar_figure
from cndpy.transform import compute_vx
from cndpy.variance import cumulative_series, get_cumulative_variance


def _pipeline(df_raw, nutrient_cols):
    df, d_plus_1 = compute_vx(df_raw, nutrient_cols)
    var_df = get_cumulative_variance(df, nutrient_cols)
    cubic_fits = fit_cubic_inflections(var_df, nutrient_cols)
    cutoff_result = select_yield_cutoff(df, cubic_fits)
    norms = compute_norms(cutoff_result['high_df'], nutrient_cols)
    prop_low = low_yield_proportion(len(df), len(cutoff_result['high_df']))
    critical_r2 = critical_chi_square(prop_low, d_plus_1)
    indices_df = compute_indices(df, nutrient_cols, norms, critical_r2)
    cum_series_by_col = {
        col: cumulative_series(var_df, f'V_{col}') for col in nutrient_cols + ['R']
    }
    return dict(
        var_df=var_df, cubic_fits=cubic_fits, cutoff_result=cutoff_result, norms=norms,
        prop_low=prop_low, critical_r2=critical_r2, d_plus_1=d_plus_1,
        indices_df=indices_df, cum_series_by_col=cum_series_by_col,
    )


def test_cumulative_variance_figure_nopal(nopal_df, nopal_nutrient_cols):
    ctx = _pipeline(nopal_df, nopal_nutrient_cols)
    best_col = select_representative_column(ctx['cubic_fits'], ctx['cutoff_result']['y_star_ref'])
    fig = cumulative_variance_figure(
        ctx['cum_series_by_col'], best_col, ctx['cubic_fits'],
        ctx['cutoff_result']['cutoff'], 'kg planta-1', nopal_nutrient_cols,
    )
    assert isinstance(fig, go.Figure)
    # 5 nutrients + R = 6 marker series, plus 1 cubic-fit overlay line
    assert len(fig.data) == 7


def test_chisquare_cdf_figure_nopal(nopal_df, nopal_nutrient_cols):
    ctx = _pipeline(nopal_df, nopal_nutrient_cols)
    fig = chisquare_cdf_figure(ctx['indices_df'], ctx['critical_r2'], ctx['prop_low'], ctx['d_plus_1'])
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 2  # theoretical line + empirical markers


def test_diagnosis_bar_figure_nopal(nopal_df, nopal_nutrient_cols):
    ctx = _pipeline(nopal_df, nopal_nutrient_cols)
    sample = {f'I_{c}': 0.5 for c in nopal_nutrient_cols + ['R']}
    fig = diagnosis_bar_figure(sample, nopal_nutrient_cols, r2_val=1.5,
                                critical_r2=ctx['critical_r2'], d_plus_1=ctx['d_plus_1'])
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 1


def test_minimal_synthetic_dataset_does_not_raise():
    df_raw = pd.DataFrame({
        'yield': [10.0, 8.0, 6.0, 4.0, 2.0, 1.0],
        'A': [1.1, 1.0, 0.9, 0.8, 0.7, 0.6],
        'B': [2.1, 2.0, 1.9, 1.8, 1.7, 1.6],
    })
    ctx = _pipeline(df_raw, ['A', 'B'])
    best_col = select_representative_column(ctx['cubic_fits'], ctx['cutoff_result']['y_star_ref'])
    fig = cumulative_variance_figure(
        ctx['cum_series_by_col'], best_col, ctx['cubic_fits'],
        ctx['cutoff_result']['cutoff'], 't ha-1', ['A', 'B'],
    )
    assert isinstance(fig, go.Figure)

    fig_chi = chisquare_cdf_figure(ctx['indices_df'], ctx['critical_r2'], ctx['prop_low'], ctx['d_plus_1'])
    assert isinstance(fig_chi, go.Figure)
