# MetaboAge: Metabolomic Aging Clock

This repository contains the code and workflow for constructing the Metabolomic Aging Clock (**MetaboAge**) using **184 NMR-based metabolic features** and **Stacked Ensemble models** in the UK Biobank.

## Overview

- **Normalization**: Box-Cox transformation to address skewed distributions.
- **Missing Value Imputation**: k-Nearest Neighbors (k=10).
- **Outlier Detection**: Mahalanobis distance in PCA space.
- **Modeling**: Stacked ensemble (XGBoost, LightGBM, CatBoost) with Elastic Net (α = 0.5) Regression as a meta-learner.
- **Sex-wise Modeling**: Separate models for Men and Women, addressing sex-specific metabolic aging rates. 

## Files

- `MetaboAge_feature_selection.R`: Feature selection and preprocessing.
- `MetaboAge_imputation.R`: KNN imputation for missing values.
- `MetaboAge_outlier_detection.R`: Outlier detection and exclusion.
- `MetaboAge_stacked_model.R`: Main model training and evaluation separately for Women and Men.
- `README.md`: This file.
- `LICENSE`: License

## How to Use

1. Prepare your metabolomics data in R.
2. Run the scripts in order as listed above.
3. Follow the instructions and comments in each script for detailed steps and configuration.

## Sensitivity Analysis

1. Compared stacked ensemble vs. single models (XGBoost, LightGBM, CatBoost).
2. Stacked ensemble consistently achieved highest correlation and lowest RMSE/MAE across sexes.
3. Details and results are provided in Table S1 of the manuscript. 


## Reference

If you use this code, please cite:

> Mostafaei S, et al. (2025) "Precision Prediction of Alzheimer's Disease and Related Dementias Using Integrative Multi-Omics Aging Clocks and Genetic Data" [Manuscript].  

## License

This project is licensed under the MIT License

## Contact

For questions or contributions, please contact: • Dr. Shayan Mostafaei (shayan.mostafaei@ki.se) • Dr. Sara Hägg (sara.hagg@ki.se)
