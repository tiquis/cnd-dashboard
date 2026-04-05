# CND Dashboard: An Open-Source Interactive Tool for Automated CND Norms 🧪

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge.svg)]( https://cnd-dashboard.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

This repository contains the source code for **CND Dashboard**, an interactive web application designed to automate the complete workflow of **Compositional Nutrient Diagnosis (CND)** in agricultural crops.

## 📝 Scientific Overview

Compositional Nutrient Diagnosis (CND) is a multivariate approach for diagnosing plant nutrient imbalances that accounts for the **closed nature** of tissue compositional data. By utilizing row-centered log-ratio (**clr**) transformations, CND circumvents the constraints of sub-compositional incoherence and the "closure" problem inherent in traditional dry-matter concentration percentages.

### Key Scientific Features:
- **Automated clr Transformation:** Systematic conversion of nutrient concentrations into scale-invariant centered log-ratios.
- **Iterative Cate-Nelson Partitioning:** Optimization of the high-yield subpopulation cutoff through the cumulative variance ratio function and cubic polynomial fitting.
- **Probabilistic Diagnosis:** Real-time calculation of CND indices and the Global Imbalance Index ($I_{CND}$), validated against the $\chi^2$ distribution.
- **Numerical Fidelity:** Computational engine validated with historical datasets (e.g., *Opuntia ficus-indica*), ensuring high-precision replication of established norms.
- **Professional Reporting:** Exportation of high-resolution analytical graphs and comprehensive diagnostic reports in PDF and Excel formats.

## 🚀 Installation and Local Deployment
