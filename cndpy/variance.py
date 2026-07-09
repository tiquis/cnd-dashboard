import numpy as np
import pandas as pd


def get_cumulative_variance(df: pd.DataFrame, nutrient_cols: list[str]) -> pd.DataFrame:
    """
    Eq. [8] Khiari et al. (2001):
      fi(VX) = Var(VX of n1 obs)  /  Var(VX of n2 obs)
    where data are sorted descending by yield;
    at iteration i (= n1-1), n1 = i+1 observations are in the HIGH group
    and n2 = n - n1 are in the LOW group.
    The yield_cut assigned to iteration i is the yield of the LAST (lowest)
    observation in the HIGH group (df_sorted.iloc[i], 0-indexed), matching
    Table 1 of Khiari (where the 5th iteration is labelled with 8.21 Mg ha-1).
    """
    df_sorted = df.sort_values('yield', ascending=False).reset_index(drop=True)
    n = len(df_sorted)
    results = []
    # i is the size of the HIGH group: starts at 2 (need ≥2 for variance),
    # stops at n-2 (need ≥2 in LOW group)
    for i in range(2, n - 1):
        high = df_sorted.iloc[:i]       # n1 = i obs (highest yields)
        low  = df_sorted.iloc[i:]       # n2 = n-i obs
        for vcol in [f'V_{col}' for col in nutrient_cols] + ['V_R']:
            var_n1 = high[vcol].var(ddof=1) if len(high) > 1 else np.nan
            var_n2 = low[vcol].var(ddof=1)  if len(low)  > 1 else np.nan
            # Eq. [8] verified against Table 1: fi = Var(n2=low) / Var(n1=high)
            # (The text says Var(n1)/Var(n2) but n1 is the REMAINDER in Cate-Nelson
            #  stepping — confirmed by reproducing paper's fi(VP)=129.47 at i=1)
            if pd.isna(var_n1) or pd.isna(var_n2) or var_n1 == 0:
                ratio = 0.0
            else:
                ratio = var_n2 / var_n1
            # BUG FIX (validated against the original 2004 Excel worksheet,
            # TODO-nopal.xls): the fi/fic value computed when the HIGH group has
            # just grown to size i must be plotted against the yield of the
            # *previous* boundary observation (df_sorted.iloc[i-2]), not the
            # newly-added one (iloc[i-1]). The original spreadsheet puts the
            # first ratio (i=2, top two yields) on the row of the single
            # highest yield, confirmed by reproducing fic(VN)=5.466433 etc.
            # to 6 significant figures. Using iloc[i-1] shifts every cubic fit
            # by one observation and silently corrupts the inflection points.
            results.append({
                'yield_cut': df_sorted['yield'].iloc[i - 2],
                'nut': vcol,
                'f_i': ratio
            })
    return pd.DataFrame(results)


def cumulative_series(var_df: pd.DataFrame, vcol: str) -> pd.DataFrame:
    """Aggregate get_cumulative_variance() output for one V_x column into the
    cumulative variance ratio function FiC (Eq. [9]).

    Selects this nutrient's rows in their natural iteration order (i = 2..n-2)
    WITHOUT groupby('yield_cut') — tied yields must stay as separate rows, or
    one of the n-3 iterations is silently dropped (BUG FIX, validated against
    TODO-nopal.xls). When total f_i is <= 0 (a nutrient that never discriminates
    between groups), returns a flat 0.0 series instead of dividing by zero
    (BUG FIX 1) rather than raising or producing NaN.

    Returns a DataFrame with columns ['yield_cut', 'f_i', 'FiC'].
    """
    sub = var_df[var_df['nut'] == vcol].reset_index(drop=True).copy()
    total = sub['f_i'].sum()
    if total > 0:
        sub['FiC'] = sub['f_i'].cumsum() / total * 100
    else:
        sub['FiC'] = 0.0
    return sub
