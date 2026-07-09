import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import chi2

from .norms import chisquare_goodness_of_fit

_MARKER_STYLES = {
    'fic(VN)':  dict(symbol='diamond',      color='navy',    size=9),
    'fic(VP)':  dict(symbol='square',        color='magenta', size=9),
    'fic(VK)':  dict(symbol='triangle-up',   color='gold',    size=10),
    'fic(VCa)': dict(symbol='x',             color='teal',    size=9),
    'fic(VMg)': dict(symbol='cross',         color='black',   size=9),
    'fic(VR)':  dict(symbol='circle',        color='brown',   size=10),
}


def cumulative_variance_figure(
    cum_series_by_col: dict[str, pd.DataFrame],
    best_col: str,
    cubic_fits: dict,
    cutoff: float,
    yield_unit: str,
    nutrient_cols: list[str],
) -> go.Figure:
    """Multi-nutrient plot of the cumulative variance ratio function, style
    of Khiari et al. (2001): one marker series per nutrient expression, a
    cubic-fit overlay line (with equation annotation) on the representative
    `best_col` expression, and a vertical dashed cutoff line.
    """
    nut_labels = {col: f'fic(V{col})' for col in nutrient_cols + ['R']}

    fig = go.Figure()

    for col in nutrient_cols + ['R']:
        label = nut_labels[col]
        ms = _MARKER_STYLES.get(label, dict(symbol='circle', color='gray', size=9))
        sub = cum_series_by_col.get(col)
        if sub is None or sub.empty:
            continue
        sub = sub.sort_values('yield_cut', ascending=False)
        fig.add_trace(go.Scatter(
            x=sub['yield_cut'],
            y=sub['FiC'],
            mode='markers',
            name=label,
            marker=dict(
                symbol=ms['symbol'],
                color=ms['color'],
                size=ms['size'],
                line=dict(width=1, color=ms['color'])
            )
        ))

    poly_title_line = ""
    if best_col in cubic_fits:
        coeffs_b, r2_poly, xk, _yk = cubic_fits[best_col]
        sort_idx = np.argsort(xk)[::-1]
        xk_s = xk[sort_idx]
        poly_y = np.polyval(coeffs_b, xk_s)

        a_b, b_b, c_b, d_b = coeffs_b

        def sgn(v):
            return "+" if v >= 0 else "-"

        eq_str = (
            f"F<sup>C</sup><sub>i</sub>(V<sub>{best_col}</sub>) = "
            f"{a_b:.3f}Y<sup>3</sup> "
            f"{sgn(b_b)} {abs(b_b):.2f}Y<sup>2</sup> "
            f"{sgn(c_b)} {abs(c_b):.2f}Y "
            f"{sgn(d_b)} {abs(d_b):.0f}"
        )
        poly_title_line = f"<br><sub>{eq_str}    R² = {r2_poly:.2f}</sub>"

        fig.add_trace(go.Scatter(
            x=xk_s, y=poly_y,
            mode='lines',
            name=f'Cubic fit (V{best_col})',
            line=dict(color='black', width=3),
            showlegend=True
        ))

    fig.add_vline(
        x=cutoff,
        line_dash="dash",
        line_color="red",
        line_width=2,
        annotation_text=f"  Yield cutoff = {cutoff:.2f} {yield_unit}",
        annotation_position="top left",
        annotation_font=dict(size=13)
    )

    fig.update_layout(
        title=dict(
            text=(
                f"Equations relating yield ({yield_unit}) to the<br>"
                f"cumulative variance ratio function in S<sup>{len(nutrient_cols)}</sup>"
                + poly_title_line
            ),
            x=0.5,
            xanchor='center',
            font=dict(size=18)
        ),
        height=750,
        width=1100,
        margin=dict(l=90, r=200, t=130, b=90),
        legend=dict(
            title=None,
            orientation="v",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            font=dict(size=14),
            bgcolor="rgba(255,255,255,0.95)",
            bordercolor="black",
            borderwidth=1
        ),
        xaxis=dict(
            title=dict(text=f"Yield ({yield_unit})", font=dict(size=16)),
            tickfont=dict(size=14),
            gridcolor="lightgray",
            autorange='reversed',
            showline=True, linecolor='black', mirror=True,
            zeroline=False
        ),
        yaxis=dict(
            title=dict(text="Cumulative variance ratio function (%)", font=dict(size=16)),
            tickfont=dict(size=14),
            gridcolor="lightgray",
            range=[0, 101],
            showline=True, linecolor='black', mirror=True,
            zeroline=False
        ),
        plot_bgcolor="white",
        paper_bgcolor="white"
    )

    return fig


def chisquare_cdf_figure(
    indices_df: pd.DataFrame, critical_r2: float, prop_low: float, d_plus_1: int,
) -> go.Figure:
    """Empirical vs. theoretical chi-square CDF of CND r^2, style of
    Magallanes-Quintanar et al. (2004), with goodness-of-fit R^2 computed via
    linear regression of the empirical CDF on the theoretical CDF values.
    """
    r2_obs = np.sort(indices_df['CND_r2'].dropna().values)
    n_obs_r2 = len(r2_obs)
    emp_cdf = np.arange(1, n_obs_r2 + 1) / n_obs_r2

    x_theory = np.linspace(0, max(r2_obs) * 1.15, 500)
    y_theory = chi2.cdf(x_theory, df=d_plus_1)

    r2_fit = chisquare_goodness_of_fit(indices_df, d_plus_1)

    fig_chi = go.Figure()

    fig_chi.add_trace(go.Scatter(
        x=x_theory, y=y_theory,
        mode='lines',
        name=f'χ²({d_plus_1}) CDF (theoretical)',
        line=dict(color='black', width=2)
    ))

    fig_chi.add_trace(go.Scatter(
        x=r2_obs, y=emp_cdf,
        mode='markers',
        name='CND r² (empirical CDF)',
        marker=dict(symbol='circle', color='black', size=7,
                    line=dict(width=1, color='black'))
    ))

    fig_chi.add_vline(
        x=critical_r2, line_dash="dash", line_color="black", line_width=1.5
    )
    fig_chi.add_hline(
        y=1 - prop_low, line_dash="dash", line_color="black", line_width=1.5
    )

    fig_chi.add_annotation(
        x=critical_r2, y=0,
        text=f"<b>{critical_r2:.2f}</b>",
        showarrow=False, xanchor='center', yanchor='top',
        yshift=-12, font=dict(size=12)
    )
    fig_chi.add_annotation(
        x=0, y=1 - prop_low,
        text=f"<b>{1 - prop_low:.3f}</b>",
        showarrow=False, xanchor='right', yanchor='middle',
        xshift=-6, font=dict(size=12)
    )

    fig_chi.update_layout(
        title=dict(
            text=(
                f"The χ² cumulative distribution function with {d_plus_1} df<br>"
                f"<sub>R² (empirical vs. theoretical) = {r2_fit:.4f}"
                f"   •   Critical r² = {critical_r2:.4f}"
                f"   •   df = {d_plus_1}"
                f"   •   Low-yield proportion = {prop_low*100:.1f}%</sub>"
            ),
            x=0.5, xanchor='center', font=dict(size=16)
        ),
        height=620, width=800,
        margin=dict(l=80, r=60, t=120, b=80),
        xaxis=dict(
            title=dict(text="Chi-square or CND r²", font=dict(size=15)),
            tickfont=dict(size=13), gridcolor="lightgray",
            range=[0, max(r2_obs) * 1.15],
            showline=True, linecolor='black', mirror=True, zeroline=False
        ),
        yaxis=dict(
            title=dict(text="Cumulative distribution function", font=dict(size=15)),
            tickfont=dict(size=13), gridcolor="lightgray",
            range=[0, 1.05],
            showline=True, linecolor='black', mirror=True, zeroline=False
        ),
        legend=dict(
            orientation="v", yanchor="bottom", y=0.05,
            xanchor="right", x=0.98, font=dict(size=13),
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="black", borderwidth=1
        ),
        plot_bgcolor="white", paper_bgcolor="white"
    )

    return fig_chi


def diagnosis_bar_figure(
    indices: dict[str, float], nutrient_cols: list[str], r2_val: float,
    critical_r2: float, d_plus_1: int,
) -> go.Figure:
    """Bar chart of I_x for a single diagnosed sample, with the most-limiting
    nutrient (highest |I_x|) highlighted in red."""
    index_labels = [f"I({col})" for col in nutrient_cols + ['R']]
    index_values = [indices[f'I_{col}'] for col in nutrient_cols + ['R']]

    max_abs = max(abs(x) for x in index_values)
    bar_colors = [
        'crimson' if abs(v) == max_abs
        else ('steelblue' if v >= 0 else 'tomato')
        for v in index_values
    ]

    fig_bar = go.Figure(go.Bar(
        x=index_labels,
        y=index_values,
        marker_color=bar_colors,
        text=[f"{v:+.3f}" for v in index_values],
        textposition='outside',
    ))
    fig_bar.add_hline(y=0, line_color='black', line_width=1)
    fig_bar.update_layout(
        title=dict(
            text=(f"CND Nutrient Indices — r² = {r2_val:.3f}  |  "
                  f"Critical = {critical_r2:.4f}  |  df = {d_plus_1}"),
            x=0.5, xanchor='center', font=dict(size=14)
        ),
        xaxis=dict(title="Nutrient index", tickfont=dict(size=13)),
        yaxis=dict(
            title="Index value (I_X)",
            tickfont=dict(size=13),
            zeroline=False,
            gridcolor="lightgray",
        ),
        plot_bgcolor="white", paper_bgcolor="white",
        height=420,
        showlegend=False,
        margin=dict(t=60, b=60, l=60, r=40),
    )
    return fig_bar
