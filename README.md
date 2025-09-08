# MetaboAge: Metabolomic Aging Clock

This repository contains the code and workflow for constructing the Metabolomic Aging Clock (MetaboAge) using 184 NMR-based metabolic features and Stacked ensemble models in the UK Biobank.

## Overview

- **Normalization**: Box-Cox transformation to address skewed distributions.
- **Missing Value Imputation**: k-Nearest Neighbors (k=10).
- **Outlier Detection**: Mahalanobis distance in PCA space.
- **Modeling**: Stacked ensemble (XGBoost, Random Forest, Decision Tree) with Elastic Net Regression as a meta-learner.
- **Subgroup Modeling**: Separate models for age (Under 50, 50-59, 60 and older) and sex groups. 

## Files

- `MetaboAge_feature_selection.R`: Feature selection and preprocessing.
- `MetaboAge_imputation.R`: KNN imputation for missing values.
- `MetaboAge_outlier_detection.R`: Outlier detection and exclusion.
- `MetaboAge_stacked_model.R`: Main model training and evaluation.
- `MetaboAge_groupwise_models.R`: Age/sex subgroup modeling.
- `README.md`: This file.
- `LICENSE`: License

## How to Use

1. Prepare your metabolomics data in R.
2. Run the scripts in order as listed above.
3. Follow instructions in each script for details.

## Reference

If you use this code, please cite:

> Mostafaei S, et al. (2025) "Precision Prediction of Alzheimer's Disease and Related Dementias Using Integrative Multi-Omics Aging Clocks and Genetic Data" [Manuscript].  

## License

This project is licensed under the MIT License

## Contact

For questions or contributions, please contact: • Dr. Shayan Mostafaei (shayan.mostafaei@ki.se) • Dr. Sara Hägg (sara.hagg@ki.se)
