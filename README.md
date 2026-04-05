# CND Dashboard: An Open-Source Interactive Tool for Automated CND Norms 🧪

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge.svg)](https://cnd-dashboard.streamlit.app/)
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

To execute the dashboard in a local environment:

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/your-username/cnd-dashboard.git](https://github.com/your-username/cnd-dashboard.git)
   cd cnd-dashboard
2.	**Configure the environment:**
It is recommended to use a virtual environment (venv or conda).
Bash
pip install -r requirements.txt

Launch the application:
Bash
streamlit run CND_V14.py


📊 Methodology & Validation
The algorithm implements the foundational protocols of Aitchison (1986) for compositional data and the diagnostic frameworks of Parent and Lekshmanan (2010). A significant contribution of this tool is the automated detection of the "plateau" in the variance ratio function, which traditionally required manual graphical interpretation and was prone to researcher subjectivity.
The dashboard's accuracy has been verified by replicating the norms for Opuntia ficus-indica (nopal) with five-decimal precision, correcting previous inconsistencies found in legacy literature regarding subpopulation proportions (e.g., the 12/40 vs 70% discrepancy).

📄 Associated Publication
If you use this tool in your research, please cite the following work:
Magallanes-Quintanar, R., & Valdez-Cepeda, R. D. (2026). CND Dashboard: An Open-Source Interactive Tool for Automated Computation of Compositional Nutrient Diagnosis Norms and Real-Time Foliar Diagnosis in Crops. Submitted to Computers and Electronics in Agriculture.

🛠️ Technical Stack
•	Backend: Python 3.x
•	Frontend: Streamlit Framework
•	Data Orchestration: Pandas & NumPy
•	Statistical Engine: SciPy Stats
•	Visualization: Plotly Interactive Library
•	Document Generation: ReportLab (PDF Engine)

🤝 Contributions and Feedback
Academic collaboration is highly encouraged. If you wish to implement new log-ratio transformations (e.g., ilr) or provide specific crop norms for the built-in database, please open an Issue or submit a Pull Request.

⚖️ License
This project is licensed under the MIT License:
Copyright (c) 2026 Rafael Magallanes-Quintanar
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
________________________________________
Corresponding Author: tiquis@uaz.edu.mx
Institutions: Universidad Autónoma de Zacatecas (UAZ) | Universidad Autónoma Chapingo (UACh)
