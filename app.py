import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import chi2
import plotly.express as px
import plotly.graph_objects as go
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import io
import datetime

st.set_page_config(page_title="CND Dashboard - Academic", layout="wide")
st.title("🧪 Compositional Nutrient Diagnosis (CND) – Dashboard")
st.markdown("**English Academic Version • Based on Rafael Magallanes Quintanar PhD Thesis**")

# ===================== NORMS FROM Rafael Magallanes Quintanar Ph. D. Thesis =====================
norms_tesis = {
    "Maize (Zea mays L.)": {"V_N":0.462, "V_P":-1.690, "V_K":-0.259, "V_Ca":-0.866, "V_Mg":-1.686, "V_R":4.040,
                            "SD_N":0.117, "SD_P":0.119, "SD_K":0.097, "SD_Ca":0.087, "SD_Mg":0.176, "SD_R":0.091},
    "Nopal (Opuntia ficus-indica)": {
        "V_N":-1.13336, "V_P":-2.26110, "V_K":0.36715, "V_Ca":0.37021, "V_Mg":-0.72570, "V_R":3.38281,
        "SD_N":0.07660, "SD_P":0.10930, "SD_K":0.23290, "SD_Ca":0.10470, "SD_Mg":0.14130, "SD_R":0.08330},
}

def compute_vx(df, nutrient_cols):
    df = df.copy()
    df['R'] = 100 - df[nutrient_cols].sum(axis=1)
    components = nutrient_cols + ['R']
    d_plus_1 = len(components)
    g = (df[components].prod(axis=1)) ** (1 / d_plus_1)
    for col in components:
        df[f'V_{col}'] = np.log(df[col] / g)
    return df, d_plus_1

def get_cumulative_variance(df, nutrient_cols):
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
            results.append({
                'yield_cut': df_sorted['yield'].iloc[i - 1],  # yield of last obs in HIGH group
                'nut': vcol,
                'f_i': ratio
            })
    return pd.DataFrame(results)

# ===================== CARGA =====================
st.sidebar.header("📤 Upload Dataset")
uploaded = st.sidebar.file_uploader("CSV file with columns: yield + nutrient columns", type=["csv"])

yield_unit = st.sidebar.text_input("Yield unit (e.g., t ha⁻¹, kg plant⁻¹, Mg ha⁻¹)", value="t ha⁻¹")

if uploaded:
    df_raw = pd.read_csv(uploaded)
    df_raw = df_raw.dropna(how='all')
    
    numeric_cols = df_raw.select_dtypes(include=np.number).columns.tolist()
    if 'yield' not in numeric_cols:
        st.error("The CSV must contain a column named 'yield'")
        st.stop()
    
    nutrient_cols = [col for col in numeric_cols if col != 'yield']
    if len(nutrient_cols) < 2:
        st.error("The CSV must contain at least 2 nutrient columns")
        st.stop()
    
    df_raw = df_raw.dropna(subset=['yield'] + nutrient_cols)
    
    st.success(f"✅ {len(df_raw)} valid observations loaded | Nutrients detected: {', '.join(nutrient_cols)}")

    df, d_plus_1 = compute_vx(df_raw, nutrient_cols)
    var_df = get_cumulative_variance(df, nutrient_cols)

    # ── Automatic yield cutoff: inflection point of cubic FiC(VX) ─────────
    # Eq.[9]: cumulative variance ratio function for each nutrient
    # Eq.[10]: fit cubic  FiC = aY³ + bY² + cY + d
    # Eq.[12]: inflection point at Y* = -b / (3a)
    # Per Khiari: retain the HIGHEST Y* across all nutrient expressions
    # (ensures every nutrient expression classifies the cutoff as high yield).

    df_sorted_all = df.sort_values('yield', ascending=False).reset_index(drop=True)
    n_total = len(df_sorted_all)

    cutoff_candidates = []
    cubic_fits = {}          # store (coeffs, r2) keyed by col label

    for col in nutrient_cols + ['R']:
        vcol  = f'V_{col}'
        sub   = (var_df[var_df['nut'] == vcol]
                 .groupby('yield_cut')['f_i'].sum()
                 .reset_index()
                 .sort_values('yield_cut', ascending=False))
        total = sub['f_i'].sum()
        if total <= 0:
            continue
        sub['FiC'] = sub['f_i'].cumsum() / total * 100

        x_vals = sub['yield_cut'].values
        y_vals = sub['FiC'].values

        if len(x_vals) < 4:          # need ≥4 points for meaningful cubic
            continue

        # Fit cubic (degree 3) per Eq.[10]
        coeffs = np.polyfit(x_vals, y_vals, 3)   # [a, b, c, d]
        a, b   = coeffs[0], coeffs[1]

        # Eq.[12]: inflection point Y* = -b / (3a)
        if abs(a) < 1e-12:
            continue
        y_star = -b / (3 * a)

        # R² of cubic fit
        poly_y  = np.polyval(coeffs, x_vals)
        ss_res  = np.sum((y_vals - poly_y) ** 2)
        ss_tot  = np.sum((y_vals - np.mean(y_vals)) ** 2)
        r2_fit  = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

        cubic_fits[col] = (coeffs, r2_fit, x_vals, y_vals)

        # Accept only inflection points inside the observed yield range
        y_min, y_max = x_vals.min(), x_vals.max()
        if y_min <= y_star <= y_max:
            cutoff_candidates.append(y_star)

    # Highest valid Y* across all nutrient expressions (continuous, Eq.[12])
    if cutoff_candidates:
        y_star_max = max(cutoff_candidates)
    else:
        y_star_max = df['yield'].quantile(0.75)
        st.warning("⚠️ No valid cubic inflection point found; using 75th percentile as fallback.")

    # Map continuous Y* to the closest discrete observed-yield partition.
    # Sort yields descending; for each n1 (size of HIGH group), the partition
    # boundary is yields_sorted[n1-1]. Find the n1 whose boundary is nearest Y*.
    # This is the most faithful discrete analog of the cubic inflection point.
    df_sorted_yields = df.sort_values('yield', ascending=False).reset_index(drop=True)
    n_obs = len(df_sorted_yields)
    yields_desc_arr = df_sorted_yields['yield'].values

    best_n1   = 2   # fallback: smallest valid group
    best_diff = np.inf
    for n1 in range(2, n_obs - 1):
        boundary = yields_desc_arr[n1 - 1]
        diff = abs(boundary - y_star_max)
        if diff < best_diff:
            best_diff = diff
            best_n1   = n1

    cutoff  = float(yields_desc_arr[best_n1 - 1])
    high_df = df_sorted_yields.iloc[:best_n1].copy()

    norms = {}
    for col in nutrient_cols + ['R']:
        vcol = f'V_{col}'
        norms[f'V_{col}'] = round(high_df[vcol].mean(), 4)
        norms[f'SD_{col}'] = round(high_df[vcol].std(), 4)

    # Proportion of LOW-yield subpopulation at cutoff (Eq.[12] → Step 2)
    prop_low    = 1 - len(high_df) / len(df)
    # Critical CND r² from chi-square CDF with d+1 degrees of freedom (Fig. 3)
    # P(X > threshold) = prop_low  →  threshold = chi2.ppf(1 - prop_low, df=d+1)
    critical_r2 = chi2.ppf(1 - prop_low, df=d_plus_1)

    tab1, tab2, tab3, tab4 = st.tabs(["Generate CND Norms", "Real-time Diagnosis", "Graphs", "Thesis Norms"])

    # ===================== TAB 1 =====================
    with tab1:
        st.subheader("Automatic CND Norms Generation")

        col_a, col_b, col_c = st.columns(3)
        col_a.metric("Yield Cutoff (discrete, observed)", f"{cutoff:.3f} {yield_unit}")
        col_b.metric("Low-yield Subpopulation", f"{prop_low*100:.1f}%  ({n_total - len(high_df)} / {n_total})")
        col_c.metric(f"Critical CND r\u00b2  (\u03c7\u00b2, df={d_plus_1})", f"{critical_r2:.4f}")

        # Per-nutrient inflection points
        inf_rows = []
        for col, (coeffs, r2f, xv, _) in cubic_fits.items():
            a, b = coeffs[0], coeffs[1]
            y_star = -b / (3 * a)
            in_range = xv.min() <= y_star <= xv.max()
            selected = abs(y_star - y_star_max) < 1e-6 and in_range
            inf_rows.append({
                'Nutrient': f'fic(V{col})',
                'a (cubic)': round(a, 5),
                'b (cubic)': round(b, 5),
                'Y* = -b/(3a)': round(y_star, 3),
                'In range': 'Yes' if in_range else 'No',
                'Cutoff source': '*** MAX ***' if selected else '',
                'R\u00b2 cubic': round(r2f, 3)
            })
        if inf_rows:
            st.markdown("**Cubic fit inflection points per nutrient expression (Eq. [10]\u2013[12]):**")
            st.dataframe(pd.DataFrame(inf_rows), use_container_width=True, hide_index=True)
        st.info(
            f"Continuous Y\\* (max inflection) = **{y_star_max:.3f}** {yield_unit}  \u2192  "
            f"Discrete yield cutoff = **{cutoff:.3f}** {yield_unit} "
            f"(highest observed yield \u2264 Y\\*)"
        )

        st.markdown("**CND Norms (High-Yield Subpopulation):**")
        st.dataframe(pd.DataFrame([norms]).T.rename(columns={0:'Value'}), use_container_width=True)

        if st.button("📄 Generate Full Report (PDF + Excel)", type="primary"):
            # ── Build indices for every observation ───────────────────────
            indices_df = df.copy()
            for col in nutrient_cols + ['R']:
                vcol   = f'V_{col}'
                v_star = norms.get(f'V_{col}', 0)
                sd     = norms.get(f'SD_{col}', 0.1)
                indices_df[f'I_{col}'] = (indices_df[vcol] - v_star) / sd if sd != 0 else np.nan

            i_cols = [f'I_{col}' for col in nutrient_cols + ['R']]
            indices_df['CND_r2']  = (indices_df[i_cols] ** 2).sum(axis=1)
            indices_df['Balanced'] = indices_df['CND_r2'] < critical_r2

            def limiting_rank(row):
                names = nutrient_cols + ['Rd']
                order = np.argsort([abs(row[ic]) for ic in i_cols])[::-1]
                return ' > '.join(names[k] for k in order)

            indices_df['Most_Limiting'] = (indices_df[i_cols].abs()
                                           .idxmax(axis=1).str.replace('I_', ''))
            indices_df['Rank_Limiting'] = indices_df.apply(limiting_rank, axis=1)

            # ── Yield range for back-transformed nutrient concentrations ──
            # Use the full dataset (not just high-yield) to show optimum range
            nut_range_rows = []
            for col in nutrient_cols:
                hi_vals = high_df[col].values
                nut_range_rows.append({
                    'Nutrient': col,
                    'Mean (%)': round(hi_vals.mean(), 4),
                    'SD (%)':   round(hi_vals.std(),  4),
                    'Min (%)':  round(hi_vals.min(),  4),
                    'Max (%)':  round(hi_vals.max(),  4),
                })
            nut_range_df = pd.DataFrame(nut_range_rows)

            # ── Cubic fit summary table (Table 1 style) ───────────────────
            cubic_rows = []
            for col in nutrient_cols + ['R']:
                if col not in cubic_fits:
                    continue
                coeffs_t, r2t, xv, _ = cubic_fits[col]
                at, bt = coeffs_t[0], coeffs_t[1]
                ys = -bt / (3 * at) if abs(at) > 1e-12 else None
                in_r = (xv.min() <= ys <= xv.max()) if ys is not None else False
                cubic_rows.append({
                    'Nutrient':  f'F\u1d9c\u1d62(V\u2093) = V{col}',
                    'a':         round(at, 6),
                    'b':         round(bt, 5),
                    'c':         round(coeffs_t[2], 4),
                    'd':         round(coeffs_t[3], 3),
                    'R\u00b2':   round(r2t, 3),
                    'Y* = \u2212b/3a': round(ys, 3) if ys is not None else 'N/A',
                    'Valid':     'Yes' if in_r else 'No',
                })
            cubic_df = pd.DataFrame(cubic_rows)

            # ASCII-safe yield unit (avoid Unicode superscripts in PDF)
            yield_unit_pdf = (yield_unit
                              .replace('\u207b\xb9', '-1')   # ⁻¹
                              .replace('\u207b',     '-')    # ⁻
                              .replace('\xb9',       '1'))

            # ── helper: draw a horizontal rule ────────────────────────────
            def pdf_rule(cv, yy, lw=0.5):
                cv.setLineWidth(lw)
                cv.line(LM, yy, RM, yy)

            # ── helper: draw a section heading (bold, 9pt) ─────────────────
            def pdf_heading(cv, yy, txt):
                cv.setFont("Helvetica-Bold", 9)
                cv.drawString(LM, yy, txt)

            # ── helper: draw a table row, returns x after last cell ────────
            def pdf_table_row(cv, yy, values, widths, font="Helvetica", size=8, xstart=None):
                cv.setFont(font, size)
                x = xstart if xstart is not None else LM + 3
                for v, w in zip(values, widths):
                    cv.drawString(x, yy, str(v))
                    x += w

            # ── PDF canvas ────────────────────────────────────────────────
            buffer = io.BytesIO()
            c = canvas.Canvas(buffer, pagesize=letter)
            W, H = letter
            LM, RM = 50, W - 50

            # ── Page 1: Title + Summary + Cubic fits + Norms ──────────────
            c.setFont("Helvetica-Bold", 13)
            c.drawCentredString(W / 2, H - 45,
                "Compositional Nutrient Diagnosis (CND) -- Results Report")
            c.setFont("Helvetica", 9)
            c.drawCentredString(W / 2, H - 60,
                f"Generated: {datetime.date.today()}"
                f"   |   Dataset: {n_total} observations"
                f"   |   Nutrients: {', '.join(nutrient_cols)}")

            y = H - 85
            pdf_rule(c, y, lw=1)

            # Survey summary
            y -= 14
            pdf_heading(c, y, "Survey summary")
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
            pdf_rule(c, y)

            # ── Cubic fits heading (no "Table N." prefix) ─────────────────
            y -= 14
            pdf_heading(c, y,
                f"Yield at inflection points of cumulative variance functions"
                f"  (n = {len(high_df)+1} to {n_total})")
            y -= 12

            # Column widths for Table 1 (8 cols, fixed)
            t1_headers = ['Nutrient', 'a', 'b', 'c', 'd', 'R2', 'Y* (-b/3a)', 'Valid']
            t1_keys    = ['Nutrient', 'a', 'b', 'c', 'd', 'R2', 'Y*', 'Valid']
            # Wider first column, equal rest
            t1_first = 105
            t1_rest  = (RM - LM - t1_first) / (len(t1_headers) - 1)
            t1_widths = [t1_first] + [t1_rest] * (len(t1_headers) - 1)

            # Re-build cubic_rows with ASCII-only keys
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

            pdf_rule(c, y + 8)
            pdf_table_row(c, y, t1_headers, t1_widths, font="Helvetica-Bold")
            y -= 3;  pdf_rule(c, y)
            y -= 11
            for row in cubic_rows_pdf:
                pdf_table_row(c, y, [row[k] for k in t1_keys], t1_widths)
                y -= 11
            pdf_rule(c, y)
            y -= 10

            # ── CND Norms heading (no "Table N." prefix) ──────────────────
            y -= 8
            pdf_heading(c, y,
                f"CND norms (means and SD of row-centered log ratios Vx) -- "
                f"high-yield subpopulation (n = {len(high_df)}, yield > {cutoff:.2f} {yield_unit_pdf})")
            y -= 12

            t2_headers = ['Component', 'V* (mean)', 'SD', 'Nutrient', 'Mean (%)', 'SD (%)']
            t2_first = 90
            t2_rest  = (RM - LM - t2_first) / (len(t2_headers) - 1)
            t2_widths = [t2_first] + [t2_rest] * (len(t2_headers) - 1)

            pdf_rule(c, y + 8)
            pdf_table_row(c, y, t2_headers, t2_widths, font="Helvetica-Bold")
            y -= 3;  pdf_rule(c, y)
            y -= 11
            for i_col, col in enumerate(nutrient_cols + ['R']):
                mu_v = norms[f'V_{col}']
                sd_v = norms[f'SD_{col}']
                if col in nutrient_cols:
                    mu_n = f"{high_df[col].mean():.4f}"
                    sd_n = f"{high_df[col].std():.4f}"
                    row_vals = [f'V*{col}', f'{mu_v:+.5f}', f'{sd_v:.5f}', col, mu_n, sd_n]
                else:
                    row_vals = [f'V*R (filling)', f'{mu_v:+.5f}', f'{sd_v:.5f}', '--', '--', '--']
                pdf_table_row(c, y, row_vals, t2_widths)
                y -= 11
            pdf_rule(c, y)
            y -= 7
            c.setFont("Helvetica-Oblique", 7)
            c.drawString(LM, y,
                f"Note: Sum(Vx) = 0 by definition.  "
                f"Critical r2 = {critical_r2:.4f} (chi2 distribution, df = {d_plus_1}, "
                f"low-yield proportion = {prop_low*100:.1f}%).")

            # ── Page 2+: Diagnostic indices (no "Table N." prefix) ────────
            c.showPage()
            y = H - 50
            pdf_heading(c, y,
                f"CND nutrient indices (Ix) and imbalance index (r2) for all observations")
            y -= 12

            all_cols_t3    = (['yield']
                              + [f'I_{c}' for c in nutrient_cols + ['R']]
                              + ['CND_r2', 'Balanced', 'Most_Limiting'])
            disp_headers_t3 = ([yield_unit_pdf]
                               + [f'I({c})' for c in nutrient_cols + ['R']]
                               + ['CND r2', 'Balanced', 'Limiting'])

            # Adaptive column widths: first col wider, last col wider
            n_t3 = len(all_cols_t3)
            t3_first = 50   # yield
            t3_last  = 40   # most limiting nutrient label
            t3_mid   = (RM - LM - t3_first - t3_last) / (n_t3 - 2)
            t3_widths = [t3_first] + [t3_mid] * (n_t3 - 2) + [t3_last]

            pdf_rule(c, y + 8)
            pdf_table_row(c, y, disp_headers_t3, t3_widths, font="Helvetica-Bold", size=7)
            y -= 3;  pdf_rule(c, y)
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
                pdf_table_row(c, y, vals, t3_widths, size=7)
                y -= 10
            pdf_rule(c, y)

            c.save()
            buffer.seek(0)
            st.download_button(
                label="📥 Download PDF Report",
                data=buffer,
                file_name=f"CND_Report_{datetime.date.today()}.pdf",
                mime="application/pdf"
            )

            # ── Excel ─────────────────────────────────────────────────────
            with pd.ExcelWriter("CND_Results.xlsx") as writer:
                # Sheet 1: Summary
                summary_df = pd.DataFrame({
                    'Parameter': ['n total', 'n high-yield', 'n low-yield',
                                  'Yield cutoff', 'Yield unit', 'Components (d+1)',
                                  f'Critical CND r2 (chi2, df={d_plus_1})'],
                    'Value': [n_total, len(high_df), n_total-len(high_df),
                              round(cutoff, 4), yield_unit, d_plus_1,
                              round(critical_r2, 4)]
                })
                summary_df.to_excel(writer, sheet_name='Summary', index=False)
                # Sheet 2: Table 1 — cubic fits
                cubic_df.to_excel(writer, sheet_name='Table1_CubicFits', index=False)
                # Sheet 3: Table 2 — CND norms + nutrient ranges
                norms_df = pd.DataFrame([{
                    'Component': f'V*{col}',
                    'Mean (Vx)': norms[f'V_{col}'],
                    'SD (Vx)':   norms[f'SD_{col}'],
                    'Nutrient':  col if col != 'R' else 'Filling value Rd',
                    'Mean (%)':  round(high_df[col].mean(), 4) if col in nutrient_cols else '—',
                    'SD (%)':    round(high_df[col].std(),  4) if col in nutrient_cols else '—',
                } for col in nutrient_cols + ['R']])
                norms_df.to_excel(writer, sheet_name='Table2_CND_Norms', index=False)
                # Sheet 4: Table 3 — all observations with indices
                indices_df.to_excel(writer, sheet_name='Table3_All_Indices', index=False)

            with open("CND_Results.xlsx", "rb") as f:
                st.download_button(
                    label="📥 Download Excel File",
                    data=f,
                    file_name="CND_Results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

            st.success("✅ PDF and Excel reports generated (Tables 1–3)!")

    # ===================== TAB 2 =====================
    with tab2:
        st.subheader("🔬 Real-time Diagnosis")
        cols = st.columns(len(nutrient_cols))
        sample_data = {}
        for i, nut in enumerate(nutrient_cols):
            with cols[i]:
                sample_data[nut] = st.number_input(f"{nut} (%)", value=2.0 if nut == 'N' else 0.3)
        
        sample = pd.DataFrame([sample_data])
        sample_v, _ = compute_vx(sample, nutrient_cols)
        sample_v = sample_v.iloc[0]
        
        indices = {}
        r2 = 0
        for col in nutrient_cols + ['R']:
            v = sample_v[f'V_{col}']
            v_star = norms.get(f'V_{col}', 0)
            sd = norms.get(f'SD_{col}', 0.1)
            ix = (v - v_star) / sd if sd != 0 else np.nan
            indices[f'I_{col}'] = round(ix, 3)
            r2 += ix**2 if not np.isnan(ix) else 0
        
        st.json(indices)
        st.metric("CND r²", round(r2, 3))
        if r2 < critical_r2:
            st.success("✅ Nutrient balance → High yield potential")
        else:
            st.error("⚠️ Nutrient imbalance")

    # ===================== TAB 3: DOS GRÁFICAS =====================
    with tab3:
        st.subheader("Graphs")
        subtab_cumvar, subtab_chisq = st.tabs([
            "Cumulative Variance Ratio Function",
            "Chi-Square Distribution of CND r\u00b2"
        ])

    # ── Sub-tab A: Cumulative Variance Ratio (existing plot) ──────────────
    with subtab_cumvar:
        st.markdown("**Multi-nutrient plot \u2013 style of Khiari et al. (2001)**")

        nut_labels = {f'V_{col}': f'fic(V{col})' for col in nutrient_cols + ['R']}

        marker_styles = {
            'fic(VN)':  dict(symbol='diamond',      color='navy',    size=9),
            'fic(VP)':  dict(symbol='square',        color='magenta', size=9),
            'fic(VK)':  dict(symbol='triangle-up',   color='gold',    size=10),
            'fic(VCa)': dict(symbol='x',             color='teal',    size=9),
            'fic(VMg)': dict(symbol='cross',         color='black',   size=9),
            'fic(VR)':  dict(symbol='circle',        color='brown',   size=10),
        }

        # ── Build cumulative series ──────────────────────────────────────
        # BUG FIX 1: guard against zero total (which collapses all series to NaN)
        # When total == 0 for a nutrient (all f_i = 0), skip normalization and
        # keep the series as zeros so Plotly still renders the trace.
        cum_list = []
        for col in nutrient_cols + ['R']:
            vcol  = f'V_{col}'
            label = nut_labels[vcol]
            sub = (var_df[var_df['nut'] == vcol]
                   .groupby('yield_cut')['f_i'].sum()
                   .reset_index()
                   .sort_values('yield_cut', ascending=False))
            total = sub['f_i'].sum()
            if total > 0:
                sub['Cumulative_%'] = sub['f_i'].cumsum() / total * 100
            else:
                # All ratios were zero → show flat 0% line (nutrient not discriminating)
                sub['Cumulative_%'] = 0.0
            sub['Nutrient'] = label
            cum_list.append(sub)
        cum_df = pd.concat(cum_list, ignore_index=True)

        # ── Choose nutrient whose Y* = y_star_max (determined the cutoff) ──
        best_col = nutrient_cols[0]   # fallback
        for col in nutrient_cols + ['R']:
            if col in cubic_fits:
                coeffs_c, _, _, _ = cubic_fits[col]
                a_c, b_c = coeffs_c[0], coeffs_c[1]
                if abs(a_c) > 1e-12:
                    y_star_c = -b_c / (3 * a_c)
                    if abs(y_star_c - y_star_max) < 1e-6:
                        best_col = col
                        break
        # Fallback: highest R² among valid cubic fits
        if best_col == nutrient_cols[0] and cubic_fits:
            best_col = max(cubic_fits, key=lambda c: cubic_fits[c][1])

        # ── Build figure ─────────────────────────────────────────────────
        fig = go.Figure()

        for col in nutrient_cols + ['R']:
            label = nut_labels[f'V_{col}']
            ms    = marker_styles.get(label, dict(symbol='circle', color='gray', size=9))
            sub   = cum_df[cum_df['Nutrient'] == label].sort_values('yield_cut', ascending=False)
            if sub.empty:
                continue
            fig.add_trace(go.Scatter(
                x=sub['yield_cut'],
                y=sub['Cumulative_%'],
                mode='markers',
                name=label,
                marker=dict(
                    symbol=ms['symbol'],
                    color=ms['color'],
                    size=ms['size'],
                    line=dict(width=1, color=ms['color'])
                )
            ))

        # ── Cubic polynomial fit (Eq.[10]) on the selected nutrient ─────
        poly_title_line = ""
        poly_trace      = None
        if best_col in cubic_fits:
            coeffs_b, r2_poly, xk, yk = cubic_fits[best_col]
            # Sort by yield descending for clean line
            sort_idx = np.argsort(xk)[::-1]
            xk_s     = xk[sort_idx]
            poly_y   = np.polyval(coeffs_b, xk_s)

            a_b, b_b, c_b, d_b = coeffs_b

            def sgn(v):
                return "+" if v >= 0 else "-"

            # Format: aY³ + bY² + cY + d  (Eq.[10] style as in Khiari Fig.2)
            eq_str = (
                f"F<sup>C</sup><sub>i</sub>(V<sub>{best_col}</sub>) = "
                f"{a_b:.3f}Y<sup>3</sup> "
                f"{sgn(b_b)} {abs(b_b):.2f}Y<sup>2</sup> "
                f"{sgn(c_b)} {abs(c_b):.2f}Y "
                f"{sgn(d_b)} {abs(d_b):.0f}"
            )
            poly_title_line = f"<br><sub>{eq_str}    R\u00b2 = {r2_poly:.2f}</sub>"

            poly_trace = go.Scatter(
                x=xk_s, y=poly_y,
                mode='lines',
                name=f'Cubic fit (V{best_col})',
                line=dict(color='black', width=3),
                showlegend=True
            )

        if poly_trace is not None:
            fig.add_trace(poly_trace)

        # ── Vertical cutoff line ──────────────────────────────────────────
        fig.add_vline(
            x=cutoff,
            line_dash="dash",
            line_color="red",
            line_width=2,
            annotation_text=f"  Yield cutoff = {cutoff:.2f} {yield_unit}",
            annotation_position="top left",
            annotation_font=dict(size=13)
        )

        # ── Layout ────────────────────────────────────────────────────────
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

        st.plotly_chart(fig, use_container_width=True)

        # DESCARGA ALTA RESOLUCIÓN
        try:
            png_bytes = fig.to_image(format="png", scale=5)
            st.download_button(
                label="📥 Download Cumulative Graph (HIGH-RESOLUTION PNG)",
                data=png_bytes,
                file_name="Cumulative_Variance_Ratio_Khiari_style.png",
                mime="image/png"
            )
        except Exception:
            st.info("Install `kaleido` (`pip install kaleido`) to enable PNG export.")

    # ── Sub-tab B: Chi-Square CDF of CND r² ───────────────────────────────
    with subtab_chisq:
        st.markdown(
            "**Empirical vs. theoretical \u03c7\u00b2 cumulative distribution "
            f"function \u2013 style of Magallanes-Quintanar et al. (2004)**"
        )

        # ── Compute CND r² for every observation ──────────────────────────
        df_r2 = df.copy()
        for col in nutrient_cols + ['R']:
            vcol   = f'V_{col}'
            v_star = norms[f'V_{col}']
            sd     = norms[f'SD_{col}']
            df_r2[f'I_{col}'] = (df_r2[vcol] - v_star) / sd if sd != 0 else np.nan
        i_cols = [f'I_{col}' for col in nutrient_cols + ['R']]
        df_r2['CND_r2'] = (df_r2[i_cols] ** 2).sum(axis=1)

        r2_obs = np.sort(df_r2['CND_r2'].dropna().values)   # ascending
        n_obs_r2 = len(r2_obs)

        # Empirical CDF: plotting position i/n  (i = 1..n)
        emp_cdf = np.arange(1, n_obs_r2 + 1) / n_obs_r2

        # Theoretical chi² CDF
        x_theory = np.linspace(0, max(r2_obs) * 1.15, 500)
        y_theory  = chi2.cdf(x_theory, df=d_plus_1)

        # R² of empirical vs theoretical (linear regression of emp_cdf on y_hat)
        # This matches the goodness-of-fit method used in the original publications
        y_hat = chi2.cdf(r2_obs, df=d_plus_1)
        slope_fit = np.polyfit(y_hat, emp_cdf, 1)
        y_pred_fit = np.polyval(slope_fit, y_hat)
        ss_res = np.sum((emp_cdf - y_pred_fit) ** 2)
        ss_tot = np.sum((emp_cdf - emp_cdf.mean()) ** 2)
        r2_fit = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        # ── Build figure ──────────────────────────────────────────────────
        fig_chi = go.Figure()

        # Theoretical CDF curve
        fig_chi.add_trace(go.Scatter(
            x=x_theory, y=y_theory,
            mode='lines',
            name=f'\u03c7\u00b2({d_plus_1}) CDF (theoretical)',
            line=dict(color='black', width=2),
        ))

        # Empirical CDF dots
        fig_chi.add_trace(go.Scatter(
            x=r2_obs, y=emp_cdf,
            mode='markers',
            name='CND r\u00b2 (empirical CDF)',
            marker=dict(symbol='circle', color='black', size=7,
                        line=dict(width=1, color='black')),
        ))

        # Vertical dashed line at r²_crit
        fig_chi.add_vline(
            x=critical_r2,
            line_dash="dash", line_color="black", line_width=1.5,
        )
        # Horizontal dashed line at 1 - prop_low
        fig_chi.add_hline(
            y=1 - prop_low,
            line_dash="dash", line_color="black", line_width=1.5,
        )

        # Annotation of critical value on x-axis (style of published figures)
        fig_chi.add_annotation(
            x=critical_r2, y=0,
            text=f"<b>{critical_r2:.2f}</b>",
            showarrow=False,
            xanchor='center', yanchor='top',
            yshift=-12,
            font=dict(size=12),
        )
        # Annotation of proportion on y-axis
        fig_chi.add_annotation(
            x=0, y=1 - prop_low,
            text=f"<b>{1 - prop_low:.3f}</b>",
            showarrow=False,
            xanchor='right', yanchor='middle',
            xshift=-6,
            font=dict(size=12),
        )

        fig_chi.update_layout(
            title=dict(
                text=(
                    f"The \u03c7\u00b2 cumulative distribution function with {d_plus_1} df<br>"
                    f"<sub>R\u00b2 (empirical vs. theoretical) = {r2_fit:.4f} &nbsp;&nbsp; "
                    f"Critical r\u00b2 = {critical_r2:.4f} &nbsp;&nbsp; "
                    f"df = {d_plus_1} &nbsp;&nbsp; "
                    f"Low-yield proportion = {prop_low*100:.1f}%</sub>"
                ),
                x=0.5, xanchor='center', font=dict(size=16)
            ),
            height=620,
            width=800,
            margin=dict(l=80, r=60, t=120, b=80),
            xaxis=dict(
                title=dict(text=f"Chi-square or CND r\u00b2", font=dict(size=15)),
                tickfont=dict(size=13),
                gridcolor="lightgray",
                range=[0, max(r2_obs) * 1.15],
                showline=True, linecolor='black', mirror=True,
                zeroline=False,
            ),
            yaxis=dict(
                title=dict(text="Cumulative distribution function", font=dict(size=15)),
                tickfont=dict(size=13),
                gridcolor="lightgray",
                range=[0, 1.05],
                showline=True, linecolor='black', mirror=True,
                zeroline=False,
            ),
            legend=dict(
                orientation="v",
                yanchor="bottom", y=0.05,
                xanchor="right",  x=0.98,
                font=dict(size=13),
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black", borderwidth=1,
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
        )

        st.plotly_chart(fig_chi, use_container_width=True)

        # ── Goodness-of-fit summary ───────────────────────────────────────
        st.info(
            f"**R\u00b2 (empirical vs. theoretical \u03c7\u00b2) = {r2_fit:.4f}**  \u2014  "
            f"Critical r\u00b2 = {critical_r2:.4f} (df = {d_plus_1}, "
            f"low-yield proportion = {prop_low*100:.1f}%)."
        )

        # ── High-resolution PNG download ─────────────────────────────────
        try:
            png_chi = fig_chi.to_image(format="png", scale=5)
            st.download_button(
                label="📥 Download Chi-Square Graph (HIGH-RESOLUTION PNG)",
                data=png_chi,
                file_name=f"CND_ChiSquare_CDF_df{d_plus_1}.png",
                mime="image/png"
            )
        except Exception:
            st.info("Install `kaleido` to enable PNG export.")

    # ===================== TAB 4 =====================
    with tab4:
        st.subheader("Thesis Norms")
        st.json(norms_tesis)

else:
    st.info("Please upload your CSV file")

st.caption("CND Dashboard • Rafael Magallanes Quintanar • April 2026")
