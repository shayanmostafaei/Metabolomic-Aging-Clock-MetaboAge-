# MetaboAge: Metabolomic Aging Clock (UK Biobank NMR)

This repository contains the code and workflow to construct **MetaboAge**, a metabolomic aging clock trained on **184 NMR-based metabolic features** from the UK Biobank and built using a **sex-stratified stacked ensemble** framework.

MetaboAge outputs a *predicted biological age* from metabolomics and can be used to derive **metabolomic age acceleration** (MetaboAge − chronological age) for downstream analyses.

---

## What this repo implements (methodology)

MetaboAge is trained using the following pipeline:

1. **Normalization / Transformation**
   - **Box–Cox transformation** to reduce skewness and improve model fit.

2. **Missing value handling**
   - **k-Nearest Neighbors (KNN) imputation with k = 9**.

3. **Outlier detection**
   - **Mahalanobis distance** using a **robust covariance estimator** (outliers are excluded prior to training).
   - (Implementation may use PCA space depending on the script version; the key requirement is robust Mahalanobis-based flagging.)

4. **Modeling**
   - **Stacked ensemble** with base learners:
     - XGBoost
     - LightGBM
     - CatBoost
   - **Elastic Net regression** as the meta-learner (**alpha = 0.5**).

5. **Sex-stratified modeling**
   - Separate models are trained for **Women** and **Men** to account for sex-specific metabolic aging patterns.

---

## Leakage-safe training for downstream ADRD prediction (recommended)

If you are using MetaboAge as a predictor in ADRD (or any downstream outcome model), we recommend:

1. **Split the analytic cohort first** (e.g., 70/30 train/test at the cohort level).
2. Train MetaboAge **only within the training set**.
3. Use **5-fold cross-validation within the training set** to generate **out-of-fold (OOF) MetaboAge predictions** for each training participant.
4. Fit the final MetaboAge models on the full training set and generate predictions for the held-out test set.

---

## Expected performance (sanity checks)

On the MetaboAge development dataset, MetaboAge typically achieves:

- **Correlation with chronological age**: ~0.74  
- **RMSE**: ~5.3 years  

In sensitivity analyses, the **stacked ensemble** consistently outperforms single base learners (XGBoost, LightGBM, CatBoost) across sexes in correlation and error metrics.

(See the associated manuscript supplementary results for details.)

---

## Repository structure

### Main scripts

- `MetaboAge_feature_selection.R`  
  Feature selection and preprocessing setup (including feature list management and transformations).

- `MetaboAge_imputation.R`  
  Missing value imputation using KNN (k=9).

- `MetaboAge_outlier_detection.R`  
  Robust Mahalanobis-based outlier detection and exclusion.

- `MetaboAge_stacked_model.R`  
  Sex-stratified stacked ensemble training, evaluation, and model saving.

### Other files

- `README.md`  
  This documentation.

- `LICENSE`  
  MIT License.

---

## Input data requirements

You must provide metabolomics data in a tabular format containing:

- Unique participant identifier (recommended)
- **Chronological age** (years)
- **Sex** (coded consistently; e.g., “Female/Male” or 0/1)
- **184 NMR metabolic features** (columns)

> Note: UK Biobank data cannot be redistributed. This repository does not include any raw UKB data.

---

## Quickstart

1. Clone this repository:
   ```bash
   git clone https://github.com/shayanmostafaei/Metabolomic-Aging-Clock-MetaboAge-
   cd Metabolomic-Aging-Clock-MetaboAge- 

## Reference

If you use this code, please cite:

> Mostafaei S, et al. (2025) "Precision Prediction of Alzheimer's Disease and Related Dementias Using Integrative Multi-Omics Aging Clocks and Genetic Data" [Manuscript].  

## Contact

For questions or contributions, please contact: • Dr. Shayan Mostafaei (shayan.mostafaei@ki.se) • Dr. Sara Hägg (sara.hagg@ki.se)
