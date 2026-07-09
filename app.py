import streamlit as st
import pandas as pd
import numpy as np
import datetime

import cndpy

st.set_page_config(page_title="CND Dashboard - Academic", layout="wide")
st.title("🧪 Compositional Nutrient Diagnosis (CND) – Dashboard")
st.markdown("**English Academic Version • Based on Rafael Magallanes Quintanar PhD Thesis**")

# ===================== NORMAS DE TU TESIS =====================
norms_tesis = {
    "Maize (Zea mays L.)": {"V_N":0.462, "V_P":-1.690, "V_K":-0.259, "V_Ca":-0.866, "V_Mg":-1.686, "V_R":4.040,
                            "SD_N":0.117, "SD_P":0.119, "SD_K":0.097, "SD_Ca":0.087, "SD_Mg":0.176, "SD_R":0.091},
    "Nopal (Opuntia ficus-indica)": {
        "V_N":-1.13336, "V_P":-2.26110, "V_K":0.36715, "V_Ca":0.37021, "V_Mg":-0.72570, "V_R":3.38281,
        "SD_N":0.07660, "SD_P":0.10930, "SD_K":0.23290, "SD_Ca":0.10470, "SD_Mg":0.14130, "SD_R":0.08330},
}

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

    df, d_plus_1 = cndpy.compute_vx(df_raw, nutrient_cols)
    var_df = cndpy.get_cumulative_variance(df, nutrient_cols)
    cubic_fits = cndpy.fit_cubic_inflections(var_df, nutrient_cols)
    cutoff_result = cndpy.select_yield_cutoff(df, cubic_fits)

    cutoff = cutoff_result['cutoff']
    high_df = cutoff_result['high_df']
    y_star_max = cutoff_result['y_star_ref']
    if cutoff_result['warning']:
        st.warning(f"⚠️ {cutoff_result['warning']}")

    n_total = len(df)
    norms = cndpy.compute_norms(high_df, nutrient_cols)
    prop_low = cndpy.low_yield_proportion(n_total, len(high_df))
    critical_r2 = cndpy.critical_chi_square(prop_low, d_plus_1)

    tab1, tab2, tab3, tab4 = st.tabs(["Generate CND Norms", "Real-time Diagnosis", "Graphs", "Thesis Norms"])

    # ===================== TAB 1 =====================
    with tab1:
        st.subheader("Automatic CND Norms Generation")

        col_a, col_b, col_c = st.columns(3)
        col_a.metric("Yield Cutoff (discrete, observed)", f"{cutoff:.3f} {yield_unit}")
        col_b.metric("Low-yield Subpopulation", f"{prop_low*100:.1f}%  ({n_total - len(high_df)} / {n_total})")
        col_c.metric(f"Critical CND r²  (χ², df={d_plus_1})", f"{critical_r2:.4f}")

        # Per-nutrient inflection points (display formatting only)
        inf_rows = []
        for col, (coeffs, r2f, xv, _) in cubic_fits.items():
            a, b = coeffs[0], coeffs[1]
            y_star = -b / (3 * a)
            in_range = xv.min() <= y_star <= xv.max()
            inf_rows.append({
                'Nutrient': f'fic(V{col})',
                'a (cubic)': round(a, 5),
                'b (cubic)': round(b, 5),
                'Y* = -b/(3a)': round(y_star, 3),
                'In range (valid)': 'Yes' if in_range else 'No',
                'Used in mean': 'Yes' if in_range else 'No',
                'R² cubic': round(r2f, 3)
            })
        if inf_rows:
            st.markdown("**Cubic fit inflection points per nutrient expression (Eq. [10]–[12]):**")
            st.dataframe(pd.DataFrame(inf_rows), use_container_width=True, hide_index=True)
        st.info(
            f"Reference Y\\* (mean of valid/in-context inflection points) = "
            f"**{y_star_max:.3f}** {yield_unit}  →  "
            f"Discrete yield cutoff = **{cutoff:.3f}** {yield_unit} "
            f"(highest observed yield ≤ Y\\*)"
        )

        st.markdown("**CND Norms (High-Yield Subpopulation):**")
        st.dataframe(pd.DataFrame([norms]).T.rename(columns={0:'Value'}), use_container_width=True)

        if st.button("📄 Generate Full Report (PDF + Excel)", type="primary"):
            indices_df = cndpy.compute_indices(df, nutrient_cols, norms, critical_r2)
            indices_df['Most_Limiting'] = cndpy.most_limiting_nutrient(indices_df, nutrient_cols)
            indices_df['Rank_Limiting'] = indices_df.apply(
                lambda row: cndpy.rank_limiting_nutrients(row, nutrient_cols), axis=1
            )

            # ── Cubic fit summary table (Table 1 style, Unicode, for display/Excel)
            cubic_rows = []
            for col in nutrient_cols + ['R']:
                if col not in cubic_fits:
                    continue
                coeffs_t, r2t, xv, _ = cubic_fits[col]
                at, bt = coeffs_t[0], coeffs_t[1]
                ys = -bt / (3 * at) if abs(at) > 1e-12 else None
                in_r = (xv.min() <= ys <= xv.max()) if ys is not None else False
                cubic_rows.append({
                    'Nutrient':  f'Fᶜᵢ(Vₓ) = V{col}',
                    'a':         round(at, 6),
                    'b':         round(bt, 5),
                    'c':         round(coeffs_t[2], 4),
                    'd':         round(coeffs_t[3], 3),
                    'R²':   round(r2t, 3),
                    'Y* = −b/3a': round(ys, 3) if ys is not None else 'N/A',
                    'Valid':     'Yes' if in_r else 'No',
                })
            cubic_df = pd.DataFrame(cubic_rows)

            pdf_bytes = cndpy.build_pdf_report(
                n_total=n_total, high_df=high_df, nutrient_cols=nutrient_cols,
                cutoff=cutoff, yield_unit=yield_unit, d_plus_1=d_plus_1,
                critical_r2=critical_r2, cubic_fits=cubic_fits, norms=norms,
                indices_df=indices_df, prop_low=prop_low,
            )
            st.download_button(
                label="📥 Download PDF Report",
                data=pdf_bytes,
                file_name=f"CND_Report_{datetime.date.today()}.pdf",
                mime="application/pdf"
            )

            xlsx_bytes = cndpy.build_excel_report(
                n_total=n_total, high_df=high_df, nutrient_cols=nutrient_cols,
                cutoff=cutoff, yield_unit=yield_unit, d_plus_1=d_plus_1,
                critical_r2=critical_r2, cubic_df=cubic_df, norms=norms,
                indices_df=indices_df,
            )
            st.download_button(
                label="📥 Download Excel File",
                data=xlsx_bytes,
                file_name="CND_Results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.success("✅ PDF and Excel reports generated (Tables 1–3)!")

    # ===================== TAB 2 =====================
    with tab2:
        st.subheader("🔬 Real-time Diagnosis")
        st.markdown(
            "Enter foliar nutrient concentrations (% dry matter) for the sample to diagnose. "
            "Default values are the **dataset means** — a balanced specimen near the "
            "population centre should give a low CND r²."
        )

        dataset_means = {col: round(float(df_raw[col].mean()), 3) for col in nutrient_cols}

        cols_input = st.columns(len(nutrient_cols))
        sample_data = {}
        for i, nut in enumerate(nutrient_cols):
            with cols_input[i]:
                sample_data[nut] = st.number_input(
                    f"{nut} (%)",
                    value=dataset_means[nut],
                    min_value=0.001,
                    step=0.001,
                    format="%.3f"
                )

        sample_df = pd.DataFrame([sample_data])
        sample_vx, _ = cndpy.compute_vx(sample_df, nutrient_cols)
        sample_indices_df = cndpy.compute_indices(sample_vx, nutrient_cols, norms, critical_r2)

        r2_val = float(sample_indices_df['CND_r2'].iloc[0])
        balanced = bool(sample_indices_df['Balanced'].iloc[0])
        indices = {
            f'I_{col}': round(float(sample_indices_df[f'I_{col}'].iloc[0]), 4)
            for col in nutrient_cols + ['R']
        }

        m1, m2, m3 = st.columns(3)
        m1.metric("CND r²", f"{r2_val:.3f}")
        m2.metric("Critical r² (χ², df={})".format(d_plus_1), f"{critical_r2:.4f}")
        m3.metric("Diagnosis", "✅ Balanced" if balanced else "⚠️ Imbalanced")

        if balanced:
            st.success(f"✅ Nutrient balance — high yield potential  "
                       f"(r² = {r2_val:.3f} < {critical_r2:.4f})")
        else:
            st.error(f"⚠️ Nutrient imbalance  "
                     f"(r² = {r2_val:.3f} ≥ {critical_r2:.4f})")

        st.divider()

        fig_bar = cndpy.diagnosis_bar_figure(indices, nutrient_cols, r2_val, critical_r2, d_plus_1)
        st.plotly_chart(fig_bar, use_container_width=True)

        rank_data = sorted(
            [(f"I({col})", indices[f'I_{col}']) for col in nutrient_cols + ['R']],
            key=lambda x: abs(x[1]), reverse=True
        )
        rank_df = pd.DataFrame(rank_data, columns=['Index', 'Value'])
        rank_df.insert(0, 'Rank', range(1, len(rank_df)+1))
        rank_df['|Value|'] = rank_df['Value'].abs().round(4)
        rank_df['Value'] = rank_df['Value'].round(4)
        rank_df['Status'] = rank_df['Value'].apply(
            lambda v: 'Excess' if v > 0 else 'Deficiency'
        )

        st.markdown("**Limiting-nutrient ranking** (highest |I_X| = most limiting):")
        st.dataframe(rank_df, use_container_width=True, hide_index=True)

        with st.expander("Raw index values (JSON)"):
            st.json(indices)

    # ===================== TAB 3: DOS SUBPESTAÑAS ========================
    with tab3:
        st.subheader("Graphs")
        subtab_cumvar, subtab_chisq = st.tabs([
            "Cumulative Variance Ratio Function",
            "Chi-Square Distribution of CND r²"
        ])

    # ── Sub-tab A: Cumulative Variance Ratio ──────────────────────────────
    with subtab_cumvar:
        st.markdown("**Multi-nutrient plot – style of Khiari et al. (2001)**")

        cum_series_by_col = {
            col: cndpy.cumulative_series(var_df, f'V_{col}') for col in nutrient_cols + ['R']
        }
        best_col = cndpy.select_representative_column(cubic_fits, y_star_max)

        fig = cndpy.cumulative_variance_figure(
            cum_series_by_col, best_col, cubic_fits, cutoff, yield_unit, nutrient_cols
        )
        st.plotly_chart(fig, use_container_width=True)

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

    # ── Sub-tab B: Chi-Square CDF of CND r² ──────────────────────────────
    with subtab_chisq:
        st.markdown(
            "**Empirical vs. theoretical χ² cumulative distribution "
            "function – style of Magallanes-Quintanar et al. (2004)**"
        )

        df_r2 = cndpy.compute_indices(df, nutrient_cols, norms, critical_r2)

        fig_chi = cndpy.chisquare_cdf_figure(df_r2, critical_r2, prop_low, d_plus_1)
        st.plotly_chart(fig_chi, use_container_width=True)

        r2_fit = cndpy.chisquare_goodness_of_fit(df_r2, d_plus_1)

        st.info(
            f"**R² (empirical vs. theoretical χ²) = {r2_fit:.4f}**"
            f"  —  Critical r² = {critical_r2:.4f}"
            f"  (df = {d_plus_1}, low-yield proportion = {prop_low*100:.1f}%)."
        )

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

st.caption("CND Dashboard 2022 Rafael Magallanes Quintanar 2022 April 2026")
