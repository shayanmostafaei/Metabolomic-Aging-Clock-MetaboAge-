# ============================================================
# MetaboAge_outlier_detection.R
# Outlier detection using ROBUST Mahalanobis distance (per manuscript)
# Input  : results/metaboage_step2_imputation/metaboage_step2_imputed.rds
# Output : results/metaboage_step3_outliers/metaboage_step3_cleaned_no_outliers.rds
# Used by: MetaboAge_stacked_model.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stats)
  library(robustbase)  # covMcd() robust covariance estimator
})

# --------------------------
# USER SETTINGS 
# --------------------------

set.seed(20250101)

IN_DIR    <- "results/metaboage_step2_imputation"
INPUT_RDS <- file.path(IN_DIR, "metaboage_step2_imputed.rds")

OUT_DIR <- "results/metaboage_step3_outliers"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata columns to keep unchanged (edit to match your dataset)
ID_COL  <- "eid"
SEX_COL <- "sex"
AGE_COL <- "age"

# Outlier detection settings
USE_PCA_SPACE <- TRUE     # TRUE matches your original approach; set FALSE to use raw feature space
PCA_VAR_KEEP  <- 0.95     # keep enough PCs to explain this fraction of variance (only if USE_PCA_SPACE = TRUE)
ALPHA_CUTOFF  <- 0.95     # chi-square cutoff (0.95 = conservative)

# --------------------------
# LOAD DATA (AFTER IMPUTATION)
# --------------------------

final_data_imputed <- readRDS(INPUT_RDS)

required_cols <- c(ID_COL, SEX_COL, AGE_COL)
missing_cols <- setdiff(required_cols, names(final_data_imputed))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "),
       "\nEdit ID_COL/SEX_COL/AGE_COL to match your data.")
}

meta_cols <- required_cols
feature_cols <- setdiff(names(final_data_imputed), meta_cols)

# Use numeric features only
feature_cols <- feature_cols[sapply(final_data_imputed[, feature_cols, drop = FALSE], is.numeric)]
if (length(feature_cols) < 2) stop("Need at least 2 numeric features for Mahalanobis outlier detection.")

X <- final_data_imputed[, feature_cols, drop = FALSE]

# Safety: if any NA remain, stop (imputation step should have removed all)
if (anyNA(X)) stop("NA values detected after imputation. Fix imputation before outlier detection.")

# --------------------------
# ROBUST MAHALANOBIS DISTANCE
# --------------------------

# Robust location + covariance (Minimum Covariance Determinant)
mcd <- covMcd(Z)
center_rob <- mcd$center
cov_rob    <- mcd$cov

md <- mahalanobis(Z, center = center_rob, cov = cov_rob)

# Chi-square threshold
df_md <- ncol(Z)
threshold <- qchisq(ALPHA_CUTOFF, df = df_md)

outlier_flag <- md > threshold
outliers_idx <- which(outlier_flag)

# --------------------------
# EXCLUDE OUTLIERS
# --------------------------

cleaned_data <- final_data_imputed[!outlier_flag, , drop = FALSE]

# --------------------------
# SAVE OUTPUTS + REPORT
# --------------------------

# Report table (one row per sample)
outlier_report <- final_data_imputed %>%
  select(all_of(meta_cols)) %>%
  mutate(
    mahalanobis_distance = as.numeric(md),
    threshold = threshold,
    df = df_md,
    alpha = ALPHA_CUTOFF,
    space = space_used,
    is_outlier = outlier_flag
  )

write.csv(outlier_report, file.path(OUT_DIR, "metaboage_step3_outlier_report.csv"), row.names = FALSE)

saveRDS(cleaned_data, file.path(OUT_DIR, "metaboage_step3_cleaned_no_outliers.rds"))

# Save artifacts useful for reproducibility/auditing
artifacts <- list(
  use_pca_space = USE_PCA_SPACE,
  pca_var_keep = PCA_VAR_KEEP,
  alpha_cutoff = ALPHA_CUTOFF,
  df = df_md,
  threshold = threshold,
  robust_center = center_rob,
  robust_cov = cov_rob,
  space_used = space_used
)
saveRDS(artifacts, file.path(OUT_DIR, "metaboage_step3_outlier_artifacts.rds"))

cat("\nDONE ✅ Outlier detection completed.\n")
cat("Input :", INPUT_RDS, "\n")
cat("Space :", space_used, "\n")
cat("Alpha :", ALPHA_CUTOFF, "| df:", df_md, "| threshold:", round(threshold, 4), "\n")
cat("Outliers detected:", length(outliers_idx), "of", nrow(final_data_imputed), "\n")
cat("Output:", file.path(OUT_DIR, "metaboage_step3_cleaned_no_outliers.rds"), "\n")
cat("Report:", file.path(OUT_DIR, "metaboage_step3_outlier_report.csv"), "\n\n")
