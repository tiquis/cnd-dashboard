import os

from cndpy.cutoff import fit_cubic_inflections, select_yield_cutoff
from cndpy.diagnosis import compute_indices, most_limiting_nutrient, rank_limiting_nutrients
from cndpy.norms import compute_norms, critical_chi_square, low_yield_proportion
from cndpy.reports import build_excel_report, build_pdf_report
from cndpy.transform import compute_vx
from cndpy.variance import get_cumulative_variance


def _build_nopal_pipeline(nopal_df, nopal_nutrient_cols):
    df, d_plus_1 = compute_vx(nopal_df, nopal_nutrient_cols)
    var_df = get_cumulative_variance(df, nopal_nutrient_cols)
    cubic_fits = fit_cubic_inflections(var_df, nopal_nutrient_cols)
    cutoff_result = select_yield_cutoff(df, cubic_fits)
    norms = compute_norms(cutoff_result['high_df'], nopal_nutrient_cols)
    prop_low = low_yield_proportion(len(df), len(cutoff_result['high_df']))
    critical_r2 = critical_chi_square(prop_low, d_plus_1)
    indices_df = compute_indices(df, nopal_nutrient_cols, norms, critical_r2)
    indices_df['Most_Limiting'] = most_limiting_nutrient(indices_df, nopal_nutrient_cols)
    indices_df['Rank_Limiting'] = indices_df.apply(
        lambda row: rank_limiting_nutrients(row, nopal_nutrient_cols), axis=1
    )

    cubic_rows = []
    for col in nopal_nutrient_cols + ['R']:
        if col not in cubic_fits:
            continue
        coeffs_t, r2t, xv, _ = cubic_fits[col]
        at, bt = coeffs_t[0], coeffs_t[1]
        ys = -bt / (3 * at) if abs(at) > 1e-12 else None
        cubic_rows.append({'Nutrient': col, 'a': at, 'b': bt, 'R2': r2t, 'Y*': ys})
    import pandas as pd
    cubic_df = pd.DataFrame(cubic_rows)

    return dict(
        n_total=len(df), high_df=cutoff_result['high_df'], nutrient_cols=nopal_nutrient_cols,
        cutoff=cutoff_result['cutoff'], yield_unit='kg planta-1', d_plus_1=d_plus_1,
        critical_r2=critical_r2, cubic_fits=cubic_fits, cubic_df=cubic_df, norms=norms,
        indices_df=indices_df, prop_low=prop_low,
    )


def test_build_pdf_report_returns_valid_pdf_bytes(nopal_df, nopal_nutrient_cols, tmp_path, monkeypatch):
    ctx = _build_nopal_pipeline(nopal_df, nopal_nutrient_cols)
    monkeypatch.chdir(tmp_path)

    pdf_kwargs = {k: v for k, v in ctx.items() if k != 'cubic_df'}
    pdf_bytes = build_pdf_report(**pdf_kwargs)

    assert isinstance(pdf_bytes, bytes)
    assert len(pdf_bytes) > 0
    assert pdf_bytes.startswith(b'%PDF')
    assert os.listdir(tmp_path) == []  # no disk side effects


def test_build_excel_report_returns_valid_xlsx_bytes(nopal_df, nopal_nutrient_cols, tmp_path, monkeypatch):
    ctx = _build_nopal_pipeline(nopal_df, nopal_nutrient_cols)
    monkeypatch.chdir(tmp_path)

    excel_kwargs = {k: v for k, v in ctx.items() if k not in ('cubic_fits', 'prop_low')}
    xlsx_bytes = build_excel_report(**excel_kwargs)

    assert isinstance(xlsx_bytes, bytes)
    assert len(xlsx_bytes) > 0
    assert xlsx_bytes[:2] == b'PK'  # xlsx is a zip archive
    assert os.listdir(tmp_path) == []  # no disk side effects
