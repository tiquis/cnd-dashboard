# CND Dashboard

**An Open-Source Interactive Tool for Automated Computation of Compositional Nutrient Diagnosis Norms and Real-Time Foliar Diagnosis in Crops**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://cnd-dashboard.streamlit.app/)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Version](https://img.shields.io/badge/Version-14-orange)

---

## Overview

CND Dashboard is a fully open-source, interactive web application that automates the complete **Compositional Nutrient Diagnosis (CND)** workflow — from raw CSV upload to publication-ready diagnostic reports — in a single interface requiring no programming knowledge.

CND is a multivariate approach for diagnosing plant nutrient imbalances based on row-centred log-ratio (clr) transformations of foliar tissue composition and a probabilistic χ²-based imbalance index (Parent & Dafir, 1992; Khiari et al., 2001). Despite its theoretical advantages over simpler diagnostic systems such as CVA and DRIS, CND adoption has been limited by the computational complexity of deriving locally calibrated norms from small survey databases. This tool eliminates that barrier.

> **Live application:** [https://cnd-dashboard.streamlit.app/](https://cnd-dashboard.streamlit.app/)

---

## Features

### Tab 1 — Generate CND Norms
- Accepts any CSV file with a `yield` column and at least two nutrient columns
- Automatically computes row-centred log ratios (clr values, *V*_X) for all *d* + 1 components
- Iterates the **cumulative variance ratio function** *F*_i^C(*V*_X) following the Cate–Nelson procedure (Khiari et al., 2001)
- Fits a **cubic polynomial** to each cumulative function and locates the inflection-point yield cutoff:

$$Y^* = -\frac{b}{3a}$$

- Retains the highest valid *Y*\* across all *d* + 1 nutrient expressions
- Derives CND norms (*V*\*_X ± SD\*_X) from the high-yield subpopulation
- Computes the critical χ² threshold from the exact low-yield proportion
- Generates downloadable **PDF and Excel reports** with three publication-style tables:
  - **Table 1:** Cubic-fit coefficients (*a*, *b*, *c*, *d*), *R*², and inflection-point yield *Y*\* for all nutrient expressions
  - **Table 2:** CND norms and optimum nutrient concentration ranges for the high-yield subpopulation
  - **Table 3:** clr values, nutrient indices (*I*_X), CND *r*², balance classification, and limiting-nutrient ranking for all observations

### Tab 2 — Real-Time Diagnosis
- Manual entry of foliar nutrient concentrations (% dry matter)
- **Default values = dataset means**, providing an immediate sanity check (a specimen at the population centre should produce a low CND *r*²)
- Displays three summary metrics: CND *r*², critical *r*², and balance classification
- **Bar chart** of all *I*_X values with the most-limiting nutrient highlighted in red
- **Ranking table** of nutrients ordered by |*I*_X|, indicating excess or deficiency relative to norms

### Tab 3 — Graphs
Two sub-tabs with interactive Plotly figures and high-resolution PNG download (scale = 5):

**Sub-tab A — Cumulative Variance Ratio Function**
- Multi-nutrient scatter plot with distinct marker symbols per nutrient, following the style of Khiari et al. (2001)
- Overlaid cubic polynomial fit on the nutrient that determined *Y*_cutoff, with equation and *R*² in the title
- Reversed x-axis (highest yield on the left) and vertical dashed line at *Y*_cutoff

**Sub-tab B — Chi-Square Distribution of CND *r*²**
- Empirical CDF of all CND *r*² values overlaid on the theoretical χ²(*d*+1) CDF
- Dashed reference lines at the critical *r*² (vertical) and 1 − *p*_low (horizontal), with axis annotations
- Goodness-of-fit *R*² displayed in the figure title

### Tab 4 — Reference Norms
- Built-in published CND norms for **nopal (*Opuntia ficus-indica*)** and **maize (*Zea mays* L.)** for rapid diagnosis without uploading a local dataset

---

## CND Methodology

The implementation follows exactly the seven-step procedure described in Khiari et al. (2001) and Magallanes-Quintanar et al. (2004):

| Step | Operation |
|------|-----------|
| 1 | Filling value: *R*_d = 100 − (N + P + K + Ca + Mg) |
| 2 | Geometric mean: *G* = [N · P · K · Ca · Mg · *R*_d]^(1/6) |
| 3 | clr values: *V*_X = ln(*X* / *G*), with Σ*V*_X = 0 as internal check |
| 4 | Cumulative variance ratio *F*_i^C(*V*_X) → cubic fit → inflection point *Y*\* |
| 5 | CND norms (*V*\*_X, SD\*_X) from high-yield subpopulation |
| 6 | Nutrient indices: *I*_X = (*V*_X − *V*\*_X) / SD\*_X |
| 7 | Imbalance index: *r*² = Σ *I*²_X ~ χ²(*d*+1) |

> **Implementation note:** The variance ratio is computed as Var(*n*₂) / Var(*n*₁) — low-group over high-group variance. The expression as printed in Khiari et al. (2001) reads Var(*n*₁)/Var(*n*₂); numerical verification against Table 1 of that reference (*f*₁(*V*_P) = 129.47) confirms that the correct computational direction is the inverse.

---

## Validation

The application was validated using **36 nopal (*Opuntia ficus-indica*) observations** from a field experiment at Zacatecas, Mexico (Magallanes-Quintanar et al., 2004), comprising three nitrogen fertilization treatments measured from March 2001 to February 2002.

| Parameter | Value |
|-----------|-------|
| Dataset | 36 nopal observations, Zacatecas, Mexico |
| Nutrients | N, P, K, Ca, Mg (% dry matter) |
| Yield variable | Fresh cladode biomass (kg plant⁻¹) |
| High-yield subpopulation | 11 observations (30.6%) |
| Yield cutoff (*Y*_cutoff) | 36.5 kg plant⁻¹ |
| Low-yield proportion | 69.4% (25/36) |
| Critical χ² (*df* = 6) | 3.8687 |
| Goodness-of-fit *R*² (CDF) | 0.961 |

All six CND norms were reproduced to five significant figures:

| Nutrient | *V*\*_X | SD\*_X | Mean (%) | SD (%) |
|----------|---------|--------|----------|--------|
| N | −1.13336 | 0.07657 | 0.972 | 0.111 |
| P | −2.26110 | 0.10932 | 0.315 | 0.034 |
| K | +0.36715 | 0.23285 | 4.470 | 1.048 |
| Ca | +0.37021 | 0.10473 | 4.369 | 0.460 |
| Mg | −0.72572 | 0.14130 | 1.469 | 0.232 |
| *R*_d | +3.38281 | 0.08334 | — | — |

---

## Input Data Format

The CSV file must contain:
- A column named **`yield`** (any numeric unit; specify unit in the sidebar)
- At least **two nutrient columns** with names corresponding to the elements (e.g., `N`, `P`, `K`, `Ca`, `Mg`)

**Example (`nopal.csv`):**

```
yield,N,P,K,Ca,Mg
52.0,0.95,0.31,5.12,4.21,1.48
45.4,1.02,0.33,4.87,4.35,1.51
42.2,0.98,0.30,4.61,4.18,1.62
...
```

Missing values are dropped automatically. Column order is flexible — the tool detects nutrient columns from all numeric columns except `yield`.

---

## Installation and Local Execution

### Requirements

```
python >= 3.10
streamlit
pandas
numpy
scipy
plotly
reportlab
openpyxl
kaleido        # optional — required for high-resolution PNG export
```

### Install

```bash
git clone https://github.com/<your-username>/cnd-dashboard.git
cd cnd-dashboard
pip install -r requirements.txt
```

### Run

```bash
streamlit run CND_V14.py
```

The application opens automatically at `http://localhost:8501`.

> **PNG export:** High-resolution PNG downloads require `kaleido`. Install with `pip install kaleido`. If not installed, the download buttons are replaced by an informational message.

---

## Repository Structure

```
cnd-dashboard/
├── CND_V14.py              # Main application (Streamlit)
├── nopal.csv               # Validation dataset (36 observations)
├── requirements.txt        # Python dependencies
├── README.md               # This file
└── figures/
    ├── Cumulative_Variance_Ratio_Khiari_style.png
    └── CND_ChiSquare_CDF_df6.png
```

---

## Dependencies

| Library | Purpose |
|---------|---------|
| `streamlit` | Interactive web application framework |
| `pandas` | Data loading, cleaning, and manipulation |
| `numpy` | Numerical computation (clr, polynomial fitting) |
| `scipy` | χ² cumulative distribution function (`scipy.stats.chi2`) |
| `plotly` | Interactive figures and high-resolution static graphics |
| `reportlab` | Automated PDF report generation |
| `openpyxl` | Excel report generation |
| `kaleido` | PNG export from Plotly figures (optional) |

---

## Citation

If you use CND Dashboard in your research, please cite the accompanying paper (under review):

> Magallanes-Quintanar, R.; Valdez-Cepeda, R.D. CND Dashboard: An Open-Source Interactive Tool for Automated Computation of Compositional Nutrient Diagnosis Norms and Real-Time Foliar Diagnosis in Crops. *Computers and Electronics in Agriculture* (under review).

And the methodological references:

> Khiari, L.; Parent, L.E.; Tremblay, N. Selecting the high-yield subpopulation for diagnosing nutrient imbalance in crops. *Agronomy Journal* **2001**, *93*, 802–808. https://doi.org/10.2134/agronj2001.933802x

> Magallanes-Quintanar, R.; Valdez-Cepeda, R.D.; et al. Compositional Nutrient Diagnosis in nopal (*Opuntia ficus-indica*). *Journal of the Professional Association for Cactus Development* **2004**, *6*, 78–89.

> Parent, L.E.; Dafir, M. A theoretical concept of compositional nutrient diagnosis. *Journal of the American Society for Horticultural Science* **1992**, *117*, 239–242. https://doi.org/10.21273/JASHS.117.2.239

---

## Authors

**Rafael Magallanes-Quintanar** (corresponding author)  
Unidad Académica de Ingeniería Eléctrica, Universidad Autónoma de Zacatecas  
Jardín Juárez 147, Centro, Zacatecas, C.P. 98000, México  
✉ tiquis@uaz.edu.mx

**Ricardo David Valdez-Cepeda**  
Centro Regional Universitario Centro-Norte, Universidad Autónoma Chapingo  
Apdo. Postal 196, C.P. 98001 Zacatecas, Zac., México  
✉ vacrida@hotmail.com

---

## License

This project is released under the [MIT License](LICENSE).  
The nopal validation dataset (`nopal.csv`) is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Acknowledgements

This work is based on the PhD thesis of R. Magallanes-Quintanar (Universidad Autónoma de Zacatecas) and the CND methodology developed by L.E. Parent, L. Khiari, and N. Tremblay. The authors thank the Centro Regional Universitario Centro-Norte, Universidad Autónoma Chapingo, Zacatecas, for access to the experimental field data.
