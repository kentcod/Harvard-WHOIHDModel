# Spatial Variation in Ischemic Heart Disease Across Europe (1985–2020)

This repository contains code, data, and supporting materials for an interdisciplinary spatial-lag regression study of ischemic heart disease (IHD) mortality across 53 European and Western Asian countries over a 36-year period (1985–2020). The research explores spatial autocorrelation in IHD outcomes and identifies environmental, demographic, and healthcare-related predictors.

---

## 📌 Study Objective

To quantify spatiotemporal associations between national IHD mortality (deaths per 100,000) and environmental, socioeconomic, and healthcare factors, while accounting for spatial spillover effects between neighboring countries.

---

## 📊 Key Findings

- A spatial autoregressive parameter **λ = 0.259** (p < 0.001) was the strongest predictor in the final model.
- Higher IHD mortality was significantly associated with:
  - Lower GDP per capita (β = –0.142)
  - Fewer physicians per 100k (β = –0.084)
  - Older populations (β = 0.110)
- Surprisingly, fruit and vegetable availability had a **positive** association (β = 0.080), possibly confounded by underlying health system variables.

---

## 🛠️ Repository Contents

- `spatial_model.R`: Spatial-lag regression model using `spml()` from the **splm** package, including drop-one likelihood ratio tests for predictor selection.
- `data/`: Cleaned panel dataset (1985–2020) and spatial weights matrix (kNN, k=4).
- `results/`: Model diagnostics, residual tests, and output tables.
- `.Rproj`, `.gitignore`, etc.

---

## 🔍 Supporting Workflows

- 🐍 **Preliminary data cleaning, EDA, and permutation importance modeling in Python:**  
  [Google Drive Folder](https://drive.google.com/drive/u/0/folders/11MoZeozmnmw7EGVAgh7-9OGg0eAruyXh)

- 🌍 **Geospatial visualizations and mapping in ArcGIS Pro:**  
  [Final ArcGIS StoryMap](https://arcg.is/0KCbje)

---

## 🧪 Methods Summary

- Spatial panel model with fixed effects (conditional likelihood specification)
- Spatial weights: Queen contiguity → balanced via kNN (k=4)
- Imputation: Iterative, regionally guided imputation for WHO data; kNN for weakly correlated indicators
- Model selection via drop-one likelihood ratio tests
- Residual diagnostics using Moran’s I and heteroskedasticity plots

---

## 🧬 Software Used

- **R**: `splm`, `spdep`, `sf`, `ggplot2`, `plm`, `dplyr`, `readr`, `gt`
- **Python**: `pandas`, `scikit-learn`, `matplotlib`, `seaborn`
- **ArcGIS Pro**: For spatial weights, map design, and projections

