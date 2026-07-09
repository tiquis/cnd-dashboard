import datetime
import io

import numpy as np
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas


def _pdf_rule(cv, LM, RM, yy, lw=0.5):
    cv.setLineWidth(lw)
    cv.line(LM, yy, RM, yy)


def _pdf_heading(cv, LM, yy, txt):
    cv.setFont("Helvetica-Bold", 9)
    cv.drawString(LM, yy, txt)


def _pdf_table_row(cv, LM, yy, values, widths, font="Helvetica", size=8, xstart=None):
    cv.setFont(font, size)
    x = xstart if xstart is not None else LM + 3
    for v, w in zip(values, widths):
        cv.drawString(x, yy, str(v))
        x += w


def build_pdf_report(
    *,
    n_total: int,
    high_df: pd.DataFrame,
    nutrient_cols: list[str],
    cutoff: float,
    yield_unit: str,
    d_plus_1: int,
    critical_r2: float,
    cubic_fits: dict,
    norms: dict,
    indices_df: pd.DataFrame,
    prop_low: float,
) -> bytes:
    """Multi-page CND results PDF report (title, survey summary, Table 1
    cubic fits, Table 2 CND norms, Table 3 all-observation indices).
    Pure function: no Streamlit calls, no disk writes — returns raw PDF bytes.
    """
    # ASCII-safe yield unit (avoid Unicode superscripts in PDF)
    yield_unit_pdf = (yield_unit
                      .replace('⁻\xb9', '-1')   # ⁻¹
                      .replace('⁻',     '-')    # ⁻
                      .replace('\xb9',       '1'))

    buffer = io.BytesIO()
    c = canvas.Canvas(buffer, pagesize=letter)
    W, H = letter
    LM, RM = 50, W - 50

    # ── Page 1: Title + Summary + Cubic fits + Norms ──────────────────────
    c.setFont("Helvetica-Bold", 13)
    c.drawCentredString(W / 2, H - 45,
        "Compositional Nutrient Diagnosis (CND) -- Results Report")
    c.setFont("Helvetica", 9)
    c.drawCentredString(W / 2, H - 60,
        f"Generated: {datetime.date.today()}"
        f"   |   Dataset: {n_total} observations"
        f"   |   Nutrients: {', '.join(nutrient_cols)}")

    y = H - 85
    _pdf_rule(c, LM, RM, y, lw=1)

    y -= 14
    _pdf_heading(c, LM, y, "Survey summary")
    y -= 13
    c.setFont("Helvetica", 9)
    for item in [
        f"Total observations (n): {n_total}",
        f"High-yield subpopulation (n1): {len(high_df)}  ({100 - prop_low*100:.1f}%)",
        f"Low-yield subpopulation (n2): {n_total - len(high_df)}  ({prop_low*100:.1f}%)",
        f"Yield cutoff (Y*): {cutoff:.3f} {yield_unit_pdf}",
        f"Components (d+1): {d_plus_1}",
        f"Critical CND r2 (chi2, df={d_plus_1}): {critical_r2:.4f}",
    ]:
        c.drawString(LM + 10, y, item)
        y -= 12
    y -= 4
    _pdf_rule(c, LM, RM, y)

    y -= 14
    _pdf_heading(c, LM, y,
        f"Yield at inflection points of cumulative variance functions"
        f"  (n = {len(high_df)+1} to {n_total})")
    y -= 12

    t1_headers = ['Nutrient', 'a', 'b', 'c', 'd', 'R2', 'Y* (-b/3a)', 'Valid']
    t1_keys    = ['Nutrient', 'a', 'b', 'c', 'd', 'R2', 'Y*', 'Valid']
    t1_first = 105
    t1_rest  = (RM - LM - t1_first) / (len(t1_headers) - 1)
    t1_widths = [t1_first] + [t1_rest] * (len(t1_headers) - 1)

    cubic_rows_pdf = []
    for col in nutrient_cols + ['R']:
        if col not in cubic_fits:
            continue
        coeffs_t, r2t, xv, _ = cubic_fits[col]
        at, bt = coeffs_t[0], coeffs_t[1]
        ys = -bt / (3 * at) if abs(at) > 1e-12 else None
        in_r = (xv.min() <= ys <= xv.max()) if ys is not None else False
        cubic_rows_pdf.append({
            'Nutrient': f'FiC(V) = V{col}',
            'a':  f'{at:.6f}',
            'b':  f'{bt:.5f}',
            'c':  f'{coeffs_t[2]:.4f}',
            'd':  f'{coeffs_t[3]:.3f}',
            'R2': f'{r2t:.3f}',
            'Y*': f'{ys:.3f}' if ys is not None else 'N/A',
            'Valid': 'Yes' if in_r else 'No',
        })

    _pdf_rule(c, LM, RM, y + 8)
    _pdf_table_row(c, LM, y, t1_headers, t1_widths, font="Helvetica-Bold")
    y -= 3;  _pdf_rule(c, LM, RM, y)
    y -= 11
    for row in cubic_rows_pdf:
        _pdf_table_row(c, LM, y, [row[k] for k in t1_keys], t1_widths)
        y -= 11
    _pdf_rule(c, LM, RM, y)
    y -= 10

    y -= 8
    _pdf_heading(c, LM, y,
        f"CND norms (means and SD of row-centered log ratios Vx) -- "
        f"high-yield subpopulation (n = {len(high_df)}, yield > {cutoff:.2f} {yield_unit_pdf})")
    y -= 12

    t2_headers = ['Component', 'V* (mean)', 'SD', 'Nutrient', 'Mean (%)', 'SD (%)']
    t2_first = 90
    t2_rest  = (RM - LM - t2_first) / (len(t2_headers) - 1)
    t2_widths = [t2_first] + [t2_rest] * (len(t2_headers) - 1)

    _pdf_rule(c, LM, RM, y + 8)
    _pdf_table_row(c, LM, y, t2_headers, t2_widths, font="Helvetica-Bold")
    y -= 3;  _pdf_rule(c, LM, RM, y)
    y -= 11
    for col in nutrient_cols + ['R']:
        mu_v = norms[f'V_{col}']
        sd_v = norms[f'SD_{col}']
        if col in nutrient_cols:
            mu_n = f"{high_df[col].mean():.4f}"
            sd_n = f"{high_df[col].std():.4f}"
            row_vals = [f'V*{col}', f'{mu_v:+.5f}', f'{sd_v:.5f}', col, mu_n, sd_n]
        else:
            row_vals = ['V*R (filling)', f'{mu_v:+.5f}', f'{sd_v:.5f}', '--', '--', '--']
        _pdf_table_row(c, LM, y, row_vals, t2_widths)
        y -= 11
    _pdf_rule(c, LM, RM, y)
    y -= 7
    c.setFont("Helvetica-Oblique", 7)
    c.drawString(LM, y,
        f"Note: Sum(Vx) = 0 by definition.  "
        f"Critical r2 = {critical_r2:.4f} (chi2 distribution, df = {d_plus_1}, "
        f"low-yield proportion = {prop_low*100:.1f}%).")

    # ── Page 2+: Diagnostic indices ────────────────────────────────────────
    c.showPage()
    y = H - 50
    _pdf_heading(c, LM, y,
        "CND nutrient indices (Ix) and imbalance index (r2) for all observations")
    y -= 12

    all_cols_t3    = (['yield']
                      + [f'I_{c}' for c in nutrient_cols + ['R']]
                      + ['CND_r2', 'Balanced', 'Most_Limiting'])
    disp_headers_t3 = ([yield_unit_pdf]
                       + [f'I({c})' for c in nutrient_cols + ['R']]
                       + ['CND r2', 'Balanced', 'Limiting'])

    n_t3 = len(all_cols_t3)
    t3_first = 50
    t3_last  = 40
    t3_mid   = (RM - LM - t3_first - t3_last) / (n_t3 - 2)
    t3_widths = [t3_first] + [t3_mid] * (n_t3 - 2) + [t3_last]

    _pdf_rule(c, LM, RM, y + 8)
    _pdf_table_row(c, LM, y, disp_headers_t3, t3_widths, font="Helvetica-Bold", size=7)
    y -= 3;  _pdf_rule(c, LM, RM, y)
    y -= 10

    show_df = indices_df.sort_values('yield', ascending=False).reset_index(drop=True)
    for _, row in show_df.iterrows():
        if y < 60:
            c.showPage()
            y = H - 50
        vals = []
        for dc in all_cols_t3:
            val = row[dc]
            if dc == 'Balanced':
                vals.append('Yes' if bool(val) else 'No')
            elif isinstance(val, (float, np.floating)):
                vals.append(f'{val:.3f}')
            else:
                vals.append(str(val))
        _pdf_table_row(c, LM, y, vals, t3_widths, size=7)
        y -= 10
    _pdf_rule(c, LM, RM, y)

    c.save()
    return buffer.getvalue()


def build_excel_report(
    *,
    n_total: int,
    high_df: pd.DataFrame,
    nutrient_cols: list[str],
    cutoff: float,
    yield_unit: str,
    d_plus_1: int,
    critical_r2: float,
    cubic_df: pd.DataFrame,
    norms: dict,
    indices_df: pd.DataFrame,
) -> bytes:
    """4-sheet Excel workbook (Summary, Table1_CubicFits, Table2_CND_Norms,
    Table3_All_Indices). Pure function: writes to an in-memory buffer only,
    no filesystem side effects. Returns the xlsx bytes.
    """
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        summary_df = pd.DataFrame({
            'Parameter': ['n total', 'n high-yield', 'n low-yield',
                          'Yield cutoff', 'Yield unit', 'Components (d+1)',
                          f'Critical CND r2 (chi2, df={d_plus_1})'],
            'Value': [n_total, len(high_df), n_total - len(high_df),
                      round(cutoff, 4), yield_unit, d_plus_1,
                      round(critical_r2, 4)]
        })
        summary_df.to_excel(writer, sheet_name='Summary', index=False)

        cubic_df.to_excel(writer, sheet_name='Table1_CubicFits', index=False)

        norms_df = pd.DataFrame([{
            'Component': f'V*{col}',
            'Mean (Vx)': norms[f'V_{col}'],
            'SD (Vx)':   norms[f'SD_{col}'],
            'Nutrient':  col if col != 'R' else 'Filling value Rd',
            'Mean (%)':  round(high_df[col].mean(), 4) if col in nutrient_cols else '—',
            'SD (%)':    round(high_df[col].std(),  4) if col in nutrient_cols else '—',
        } for col in nutrient_cols + ['R']])
        norms_df.to_excel(writer, sheet_name='Table2_CND_Norms', index=False)

        indices_df.to_excel(writer, sheet_name='Table3_All_Indices', index=False)

    return buffer.getvalue()
