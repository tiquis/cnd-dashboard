import numpy as np
import pandas as pd


def compute_vx(df: pd.DataFrame, nutrient_cols: list[str]) -> tuple[pd.DataFrame, int]:
    """Row-centered log-ratio (clr) transform (Aitchison, 1986).

    Adds a filling value column 'R' (Rd = 100 - sum(nutrient_cols)), the
    geometric mean g of all d+1 components per row, and one V_{col} = ln(col/g)
    column per component (including 'R'). Returns (df_with_Vx_columns, d_plus_1).
    """
    df = df.copy()
    df['R'] = 100 - df[nutrient_cols].sum(axis=1)
    components = nutrient_cols + ['R']
    d_plus_1 = len(components)
    g = (df[components].prod(axis=1)) ** (1 / d_plus_1)
    for col in components:
        df[f'V_{col}'] = np.log(df[col] / g)
    return df, d_plus_1
