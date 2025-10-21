#######################################################################################
## MetaboAge: Sex-Stratified XGBoost / LightGBM / CatBoost Stacked Ensemble Model
## Author: Shayan Mostafaei  
## GitHub: https://github.com/shayanmostafaei/Metabolomic-Aging-Clock-MetaboAge- 
#######################################################################################
# ---- 0) Setup ----
required_pkgs <- c("dplyr","impute","caret","xgboost","lightgbm","catboost",
                   "glmnet","Metrics","haven","ggplot2","tidyr","tibble")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if(pkg=="catboost") {
      message("CatBoost requires manual installation: https://catboost.ai/docs/installation/r-installation.html")
    } else install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

set.seed(123)

# ---- 1) Data Preparation ----
if(!exists("data") || !"sex" %in% names(data) || !"age" %in% names(data)) stop("Input 'data' missing required columns.")
data_ml <- data %>% filter(!is.na(age)) %>% dplyr::mutate(sex = haven::as_factor(sex) %>% droplevels())
feature_cols <- setdiff(names(data_ml), c("f.eid","sex","age"))
safe_names <- make.names(feature_cols)
names(data_ml)[match(feature_cols, names(data_ml))] <- safe_names

# Impute numeric metabolite features using KNN
data_numeric <- as.data.frame(lapply(data_ml[, safe_names], as.numeric))
imputed <- impute::impute.knn(as.matrix(data_numeric), k = 10)
data_imputed <- as.data.frame(imputed$data)
names(data_imputed) <- safe_names
data_imputed$f.eid <- as.character(data_ml$f.eid)
data_imputed$sex   <- data_ml$sex
data_imputed$age   <- data_ml$age

# ---- 2) Metrics function ----
calc_metrics <- function(truth, pred){
  valid_idx <- !is.na(pred) & !is.na(truth)
  truth <- truth[valid_idx]; pred <- pred[valid_idx]
  if(length(truth)<2) return(data.frame(R_Squared=NA, Correlation=NA, RMSE=NA, MAE=NA))
  ss_res <- sum((truth-pred)^2); ss_tot <- sum((truth-mean(truth))^2)
  data.frame(
    R_Squared = 1 - ss_res/ss_tot,
    Correlation = cor(truth, pred),
    RMSE = Metrics::rmse(truth, pred),
    MAE  = Metrics::mae(truth, pred)
  )
}

# ---- 3) Base model functions ----
fit_predict_xgb <- function(train_df, test_df, full_df=NULL, nrounds=500, eta=0.03, max_depth=6){
  feats <- setdiff(names(train_df), c("age","sex","f.eid"))
  dtrain <- xgboost::xgb.DMatrix(as.matrix(train_df[, feats]), label=train_df$age)
  dtest  <- xgboost::xgb.DMatrix(as.matrix(test_df[, feats]), label=test_df$age)
  model <- xgboost::xgb.train(params=list(objective="reg:squarederror", eta=eta, max_depth=max_depth, eval_metric="rmse"),
                              data=dtrain, nrounds=nrounds, verbose=0)
  list(model=model, pred_test=predict(model,dtest), pred_full=if(!is.null(full_df)) predict(model, xgboost::xgb.DMatrix(as.matrix(full_df[, feats]))) else NULL)
}

fit_predict_lgb <- function(train_df, test_df, full_df=NULL, nrounds=1000, learning_rate=0.03, num_leaves=31){
  feats <- setdiff(names(train_df), c("age","sex","f.eid"))
  dtrain <- lightgbm::lgb.Dataset(as.matrix(train_df[, feats]), label=train_df$age)
  val <- list(valid=lightgbm::lgb.Dataset(as.matrix(test_df[, feats]), label=test_df$age))
  params <- list(objective="regression", metric="rmse", learning_rate=learning_rate, num_leaves=num_leaves)
  model <- lightgbm::lgb.train(params, dtrain, nrounds=nrounds, valids=val, early_stopping_rounds=25, verbose=-1)
  list(model=model, pred_test=predict(model, as.matrix(test_df[, feats])), pred_full=if(!is.null(full_df)) predict(model, as.matrix(full_df[, feats])) else NULL)
}

fit_predict_cat <- function(train_df, test_df, full_df=NULL, iterations=500, learning_rate=0.03, depth=6){
  if(!requireNamespace("catboost", quietly=TRUE)){ warning("CatBoost not available"); return(list(model=NULL, pred_test=rep(NA,nrow(test_df)), pred_full=if(!is.null(full_df)) rep(NA,nrow(full_df)) else NULL)) }
  feats <- setdiff(names(train_df), c("age","sex","f.eid"))
  pool_train <- catboost::catboost.load_pool(as.matrix(train_df[, feats]), label=train_df$age)
  pool_test  <- catboost::catboost.load_pool(as.matrix(test_df[, feats]), label=test_df$age)
  params <- list(loss_function="RMSE", iterations=iterations, learning_rate=learning_rate, depth=depth, od_type="Iter", od_wait=25, verbose=FALSE)
  model <- catboost::catboost.train(pool_train, pool_test, params=params)
  list(model=model, pred_test=catboost::catboost.predict(model,pool_test), pred_full=if(!is.null(full_df)) catboost::catboost.predict(model, catboost::catboost.load_pool(as.matrix(full_df[, feats]))) else NULL)
}

# ---- 4) Holdout Split (70/30) ----
train_index <- sample(seq_len(nrow(data_imputed)), size = floor(0.7 * nrow(data_imputed)))
train_data <- data_imputed[train_index, ]
test_data  <- data_imputed[-train_index, ]

# ---- 5) Run sex-stratified base models ----
sex_levels <- levels(data_imputed$sex)
metaboage_all_predictions <- data_imputed %>% dplyr::select(f.eid, sex, age) %>%
  dplyr::mutate(MetaboAge_XGB=NA_real_, MetaboAge_LGB=NA_real_, MetaboAge_CatBoost=NA_real_, MetaboAge_Stack=NA_real_)
test_metrics_all <- list()

for(s in sex_levels){
  train_s <- train_data %>% filter(sex==s)
  test_s  <- test_data %>% filter(sex==s)
  if(nrow(train_s)<5) next
  
  xgb_res <- fit_predict_xgb(train_s, test_s, full_df=data_imputed %>% filter(sex==s))
  metaboage_all_predictions$MetaboAge_XGB[metaboage_all_predictions$sex==s] <- xgb_res$pred_full
  test_metrics_all <- append(test_metrics_all, list(calc_metrics(test_s$age, xgb_res$pred_test) %>% mutate(Sex=s, Model="XGBoost")))
  
  lgb_res <- fit_predict_lgb(train_s, test_s, full_df=data_imputed %>% filter(sex==s))
  metaboage_all_predictions$MetaboAge_LGB[metaboage_all_predictions$sex==s] <- lgb_res$pred_full
  test_metrics_all <- append(test_metrics_all, list(calc_metrics(test_s$age, lgb_res$pred_test) %>% mutate(Sex=s, Model="LightGBM")))
  
  cat_res <- fit_predict_cat(train_s, test_s, full_df=data_imputed %>% filter(sex==s))
  metaboage_all_predictions$MetaboAge_CatBoost[metaboage_all_predictions$sex==s] <- cat_res$pred_full
  test_metrics_all <- append(test_metrics_all, list(calc_metrics(test_s$age, cat_res$pred_test) %>% mutate(Sex=s, Model="CatBoost")))
}

message("Base Model Test Metrics:")
print(bind_rows(test_metrics_all))

# ---- 6) Stacked Elastic Net Meta-Learner ----
stacking_results <- list()
base_feats <- c("MetaboAge_XGB","MetaboAge_LGB","MetaboAge_CatBoost")
for(s in sex_levels){
  train_s <- train_data %>% filter(sex==s)
  test_s  <- test_data %>% filter(sex==s)
  valid_idx <- rowSums(is.na(train_s[, base_feats]))==0
  train_s <- train_s[valid_idx,]; if(nrow(train_s)<5) next
  
  x_train <- as.matrix(train_s[, base_feats]); y_train <- train_s$age
  x_test  <- as.matrix(test_s[, base_feats])
  cv <- glmnet::cv.glmnet(x_train, y_train, alpha=0.5, nfolds=5)
  meta_mod <- glmnet::glmnet(x_train, y_train, alpha=0.5, lambda=cv$lambda.min)
  
  metaboage_all_predictions$MetaboAge_Stack[metaboage_all_predictions$sex==s] <- predict(meta_mod, newx=as.matrix(metaboage_all_predictions[metaboage_all_predictions$sex==s, base_feats]), s=cv$lambda.min)
  stacking_results <- append(stacking_results, list(calc_metrics(test_s$age, as.vector(predict(meta_mod,x_test,s=cv$lambda.min))) %>% mutate(Sex=s, Model="Stacked_ElasticNet")))
}

stacking_summary <- bind_rows(stacking_results)
message("Stacked Model Test Metrics:")
print(stacking_summary)

# ---- 7) Final Evaluation and Visualization ----
chron_age <- metaboage_all_predictions$age
pred_age  <- metaboage_all_predictions$MetaboAge_Stack
plot_title <- paste0(
  "Chronological Age vs MetaboAge (Stacked Model, 70/30 Holdout)\n",
  "Correlation: ", round(cor(chron_age,pred_age,use="complete.obs"),3),
  " | R-squared: ", round(cor(chron_age,pred_age,use="complete.obs")^2,3),
  " | MAE: ", round(Metrics::mae(chron_age,pred_age),3),
  " | RMSE: ", round(Metrics::rmse(chron_age,pred_age),3)
)

p <- ggplot(metaboage_all_predictions, aes(x=age, y=MetaboAge_Stack)) +
  geom_point(size=0.5, alpha=0.8, color="#0072B2", position=position_jitter(width=0.3,height=0)) +
  geom_abline(intercept=0,slope=1,color="black",linewidth=1.2) +
  labs(title=plot_title, x="Chronological Age", y="MetaboAge (Predicted Age)") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.title=element_text(face="bold"),
        axis.text=element_text(color="black")) +
  scale_x_continuous(limits=c(38,72), breaks=seq(40,70,5)) +
  scale_y_continuous(limits=c(38,72), breaks=seq(40,70,5))
print(p)

# Save final predictions
write.csv(metaboage_all_predictions, "metaboage_all_predictions_stacked.csv", row.names=FALSE)
