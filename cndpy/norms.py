import numpy as np
import pandas as pd
from scipy.stats import chi2


def compute_norms(high_df: pd.DataFrame, nutrient_cols: list[str]) -> dict[str, float]:
    """CND norms (V*_X, SD*_X) from the high-yield subpopulation."""
    norms: dict[str, float] = {}
    for col in nutrient_cols + ['R']:
        vcol = f'V_{col}'
        norms[f'V_{col}'] = round(high_df[vcol].mean(), 5)
        norms[f'SD_{col}'] = round(high_df[vcol].std(), 5)
    return norms


def low_yield_proportion(n_total: int, n_high: int) -> float:
    """Proportion of the survey in the LOW-yield subpopulation."""
    return 1 - n_high / n_total


def critical_chi_square(prop_low: float, d_plus_1: int) -> float:
    """Critical CND r^2 from the chi-square CDF with d+1 degrees of freedom
    (Fig. 3): P(X > threshold) = prop_low -> threshold = chi2.ppf(1-prop_low, df)."""
    return chi2.ppf(1 - prop_low, df=d_plus_1)


def chisquare_goodness_of_fit(indices_df: pd.DataFrame, d_plus_1: int) -> float:
    """Goodness of fit (R^2) of the empirical CDF of CND_r2 against the
    theoretical chi-square(d_plus_1) CDF, via linear regression of the
    empirical CDF on the theoretical CDF values at each observed CND_r2."""
    r2_obs = np.sort(indices_df['CND_r2'].dropna().values)
    emp_cdf = np.arange(1, len(r2_obs) + 1) / len(r2_obs)
    y_hat = chi2.cdf(r2_obs, df=d_plus_1)
    slope_fit = np.polyfit(y_hat, emp_cdf, 1)
    y_pred = np.polyval(slope_fit, y_hat)
    ss_res = np.sum((emp_cdf - y_pred) ** 2)
    ss_tot = np.sum((emp_cdf - emp_cdf.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
