# ============================================================
# MetaboAge_feature_selection.R
# Feature selection + Box–Cox preprocessing
# Output is used by: MetaboAge_imputation.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)   # for findCorrelation() + BoxCox preProcess()
})

# --------------------------
# USER SETTINGS 
# --------------------------

set.seed(20250101)

# Input: an .rds containing a data.frame named "clean_data"
INPUT_RDS <- "data/metaboage_input_raw.rds"

# Output folder
OUT_DIR <- "results/metaboage_step1_feature_selection"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata columns to keep unchanged (edit to match your dataset)
ID_COL  <- "eid"
SEX_COL <- "sex"
AGE_COL <- "age"

# Correlation cutoff
COR_CUTOFF <- 0.99

# --------------------------
# LOAD DATA
# --------------------------

obj <- readRDS(INPUT_RDS)

# Accept either a data.frame directly OR a list containing clean_data
if (is.data.frame(obj)) {
  clean_data <- obj
} else if (!is.null(obj$clean_data) && is.data.frame(obj$clean_data)) {
  clean_data <- obj$clean_data
} else {
  stop("INPUT_RDS must contain a data.frame or a list with $clean_data (data.frame).")
}

# --------------------------
# SPLIT METADATA vs FEATURES
# --------------------------

meta_cols <- intersect(c(ID_COL, SEX_COL, AGE_COL), names(clean_data))
meta_df   <- if (length(meta_cols) > 0) clean_data[, meta_cols, drop = FALSE] else NULL

# Use numeric columns that are NOT metadata as metabolite features
feature_cols <- setdiff(names(clean_data), meta_cols)
feature_cols <- feature_cols[sapply(clean_data[, feature_cols, drop = FALSE], is.numeric)]

if (length(feature_cols) < 2) {
  stop("Not enough numeric metabolite features found. Check your input and metadata column names.")
}

X <- clean_data[, feature_cols, drop = FALSE]

# --------------------------
# QC: DROP ALL-NA / CONSTANT FEATURES
# --------------------------

n_all <- nrow(X)

all_na <- names(X)[sapply(X, function(v) all(is.na(v)))]
X <- X %>% select(-any_of(all_na))

constant <- names(X)[sapply(X, function(v) {
  vv <- v[!is.na(v)]
  length(unique(vv)) < 2
})]
X <- X %>% select(-any_of(constant))

feature_cols2 <- names(X)

# --------------------------
# REMOVE HIGHLY CORRELATED FEATURES (|r| > 0.99)
# Stable approach: caret::findCorrelation
# --------------------------

to_remove_corr <- character(0)
if (ncol(X) >= 2) {
  cor_matrix <- cor(X, use = "pairwise.complete.obs")
  to_remove_corr <- findCorrelation(cor_matrix, cutoff = COR_CUTOFF, names = TRUE, exact = TRUE)
  if (length(to_remove_corr) > 0) {
    X <- X %>% select(-any_of(to_remove_corr))
  }
}

kept_features <- names(X)

# --------------------------
# BOX–COX TRANSFORM
# Box–Cox requires strictly positive values.
# We shift ONLY features needing it and save the per-feature offsets.
# Then we fit BoxCox on the shifted data and transform.
# --------------------------

offsets <- rep(0, length(kept_features))
names(offsets) <- kept_features

mins <- sapply(X, function(v) suppressWarnings(min(v, na.rm = TRUE)))
need_shift <- names(mins)[is.finite(mins) & mins <= 0]

if (length(need_shift) > 0) {
  for (f in need_shift) {
    off <- abs(mins[[f]]) + 1
    offsets[[f]] <- off
    X[[f]] <- X[[f]] + off
  }
}

# Fit BoxCox and transform (per-feature lambdas stored in object)
bc_obj <- preProcess(X, method = c("BoxCox"))
X_bc   <- predict(bc_obj, X)

# --------------------------
# FINAL OUTPUT DATAFRAME
# --------------------------

clean_data_step1 <- if (!is.null(meta_df)) {
  bind_cols(meta_df, X_bc)
} else {
  X_bc
}

# --------------------------
# SAVE OUTPUTS FOR NEXT SCRIPTS
# --------------------------

# 1) Data after feature selection + BoxCox
saveRDS(clean_data_step1, file.path(OUT_DIR, "metaboage_step1_boxcox.rds"))

# 2) Artifacts needed to apply same transformation later (e.g., test set)
artifacts <- list(
  meta_cols = meta_cols,
  original_feature_cols = feature_cols2,
  removed_all_na = all_na,
  removed_constant = constant,
  removed_high_corr = to_remove_corr,
  kept_features = kept_features,
  cor_cutoff = COR_CUTOFF,
  offsets = offsets,
  boxcox_object = bc_obj
)

saveRDS(artifacts, file.path(OUT_DIR, "metaboage_step1_artifacts.rds"))

# 3) Simple summary
summary_tbl <- data.frame(
  step = c("start_numeric_features", "removed_all_na", "removed_constant", "removed_high_corr", "kept_features"),
  n = c(length(feature_cols), length(all_na), length(constant), length(to_remove_corr), length(kept_features))
)
write.csv(summary_tbl, file.path(OUT_DIR, "metaboage_step1_summary.csv"), row.names = FALSE)

cat("\nDONE ✅ MetaboAge feature selection + BoxCox completed.\n")
cat("Input:", INPUT_RDS, "\n")
cat("Output data:", file.path(OUT_DIR, "metaboage_step1_boxcox.rds"), "\n")
cat("Artifacts:", file.path(OUT_DIR, "metaboage_step1_artifacts.rds"), "\n")
cat("Kept features:", length(kept_features), "\n\n")
