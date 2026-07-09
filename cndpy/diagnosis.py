import numpy as np
import pandas as pd


def compute_indices(
    df_vx: pd.DataFrame,
    nutrient_cols: list[str],
    norms: dict[str, float],
    critical_r2: float,
) -> pd.DataFrame:
    """Given a dataframe that already has V_{col} columns (output of
    compute_vx), add I_{col} = (V_{col} - V*_{col}) / SD_{col} for
    col in nutrient_cols + ['R'], plus 'CND_r2' = sum(I_col^2) and
    'Balanced' = CND_r2 < critical_r2.

    Works identically for a full dataset (many rows) or a single sample
    (1-row dataframe) — this unifies what used to be two separate,
    duplicated code paths (full-report indices and real-time diagnosis).
    """
    indices_df = df_vx.copy()
    for col in nutrient_cols + ['R']:
        vcol = f'V_{col}'
        v_star = norms.get(f'V_{col}', 0)
        sd = norms.get(f'SD_{col}', 0.1)
        indices_df[f'I_{col}'] = (indices_df[vcol] - v_star) / sd if sd != 0 else np.nan

    i_cols = [f'I_{col}' for col in nutrient_cols + ['R']]
    indices_df['CND_r2'] = (indices_df[i_cols] ** 2).sum(axis=1)
    indices_df['Balanced'] = indices_df['CND_r2'] < critical_r2
    return indices_df


def rank_limiting_nutrients(indices_row: pd.Series, nutrient_cols: list[str]) -> str:
    """' > '.join of nutrient names ordered by descending |I_x|.

    Mirrors the limiting_rank() closure in the original app: the filling
    value is labelled 'Rd' here (not 'R'), matching the original display
    convention.
    """
    names = nutrient_cols + ['Rd']
    i_cols = [f'I_{col}' for col in nutrient_cols + ['R']]
    order = np.argsort([abs(indices_row[ic]) for ic in i_cols])[::-1]
    return ' > '.join(names[k] for k in order)


def most_limiting_nutrient(indices_df: pd.DataFrame, nutrient_cols: list[str]) -> pd.Series:
    """idxmax over abs(I_ columns), with the 'I_' prefix stripped."""
    i_cols = [f'I_{col}' for col in nutrient_cols + ['R']]
    return indices_df[i_cols].abs().idxmax(axis=1).str.replace('I_', '')
