# MetaboAge Missing Value Imputation (KNN)
library(VIM)
set.seed(123)
id_col <- clean_data[,1]
data_to_impute <- clean_data[,-1]
final_data_imputed <- kNN(data_to_impute, k = 10)
final_data_imputed <- cbind(id_col, final_data_imputed)
