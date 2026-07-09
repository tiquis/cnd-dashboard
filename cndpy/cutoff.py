import numpy as np
import pandas as pd

from .variance import cumulative_series

CubicFit = tuple[np.ndarray, float, np.ndarray, np.ndarray]  # (coeffs, r2, x_vals, y_vals)


def fit_cubic_inflections(
    var_df: pd.DataFrame, nutrient_cols: list[str]
) -> dict[str, CubicFit]:
    """For each nutrient column (+ 'R'), aggregate via cumulative_series(),
    fit a cubic FiC(Y) = aY^3+bY^2+cY+d (Eq. [10], np.polyfit degree 3), and
    compute R^2 of the fit. Skips columns whose total f_i <= 0 (never
    discriminates), with < 4 aggregated points, or with |a| < 1e-12
    (no meaningful inflection point).

    Returns {nutrient_col: (coeffs, r2_fit, x_vals, y_vals)}.
    """
    cubic_fits: dict[str, CubicFit] = {}

    for col in nutrient_cols + ['R']:
        vcol = f'V_{col}'
        sub = cumulative_series(var_df, vcol)
        total = sub['f_i'].sum()
        if total <= 0:
            continue

        x_vals = sub['yield_cut'].values
        y_vals = sub['FiC'].values

        if len(x_vals) < 4:          # need ≥4 points for meaningful cubic
            continue

        coeffs = np.polyfit(x_vals, y_vals, 3)   # [a, b, c, d]
        a = coeffs[0]

        if abs(a) < 1e-12:
            continue

        poly_y = np.polyval(coeffs, x_vals)
        ss_res = np.sum((y_vals - poly_y) ** 2)
        ss_tot = np.sum((y_vals - np.mean(y_vals)) ** 2)
        r2_fit = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

        cubic_fits[col] = (coeffs, r2_fit, x_vals, y_vals)

    return cubic_fits


def select_yield_cutoff(df: pd.DataFrame, cubic_fits: dict[str, CubicFit]) -> dict:
    """Map the cubic-fit inflection points (Eq. [12]: Y* = -b/(3a)) to a
    discrete Cate-Nelson yield cutoff.

    A nutrient expression is a valid ("in-context") candidate when its
    inflection point falls inside its own observed yield range; nutrients
    whose Y* extrapolates outside the data are excluded ("out of context",
    Magallanes-Quintanar et al., 2004).

    The reference cutoff y_star_ref is the MEAN of all in-context Y*
    candidates — not the single highest one — reproducing the empirical
    method actually used in Magallanes-Quintanar et al. (2004) rather than
    the literal single-maximum rule in Khiari et al. (2001).

    Returns a dict: {
        'cutoff': float, 'best_n1': int, 'high_df': pd.DataFrame,
        'y_star_ref': float, 'used_fallback': bool, 'warning': str | None,
    }
    """
    cutoff_candidates = []
    for coeffs, _r2_fit, x_vals, _y_vals in cubic_fits.values():
        a, b = coeffs[0], coeffs[1]
        y_star = -b / (3 * a)
        y_min, y_max = x_vals.min(), x_vals.max()
        if y_min <= y_star <= y_max:
            cutoff_candidates.append(y_star)

    warning = None
    used_fallback = False
    if cutoff_candidates:
        y_star_ref = float(np.mean(cutoff_candidates))
    else:
        y_star_ref = df['yield'].quantile(0.75)
        used_fallback = True
        warning = "No valid cubic inflection point found; using 75th percentile as fallback."

    df_sorted_yields = df.sort_values('yield', ascending=False).reset_index(drop=True)
    n_obs = len(df_sorted_yields)
    yields_desc_arr = df_sorted_yields['yield'].values

    best_n1 = n_obs - 2   # fallback: largest valid group
    for n1 in range(2, n_obs - 1):
        boundary = yields_desc_arr[n1 - 1]
        if boundary <= y_star_ref:
            best_n1 = n1
            break   # first match is the one with the highest boundary <= Y*

    cutoff = float(yields_desc_arr[best_n1 - 1])
    high_df = df_sorted_yields.iloc[:best_n1].copy()

    return {
        'cutoff': cutoff,
        'best_n1': best_n1,
        'high_df': high_df,
        'y_star_ref': y_star_ref,
        'used_fallback': used_fallback,
        'warning': warning,
    }


def select_representative_column(cubic_fits: dict[str, CubicFit], y_star_ref: float) -> str | None:
    """Pick the nutrient expression to display as the representative cubic
    fit overlay: the valid (in-range) candidate whose Y* is closest to the
    mean reference cutoff y_star_ref (since y_star_ref is itself a mean
    across candidates, no single nutrient matches it exactly).

    Falls back to the highest-R² fit among all cubic_fits if none are
    in-range (should rarely trigger).
    """
    best_col = None
    best_dist = None
    for col, (coeffs, _r2_fit, x_vals, _y_vals) in cubic_fits.items():
        a, b = coeffs[0], coeffs[1]
        if abs(a) <= 1e-12:
            continue
        y_star = -b / (3 * a)
        if not (x_vals.min() <= y_star <= x_vals.max()):
            continue
        dist = abs(y_star - y_star_ref)
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_col = col

    if best_col is None and cubic_fits:
        best_col = max(cubic_fits, key=lambda c: cubic_fits[c][1])

    return best_col
