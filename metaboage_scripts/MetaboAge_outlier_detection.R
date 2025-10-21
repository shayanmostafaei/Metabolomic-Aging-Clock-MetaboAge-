# MetaboAge Outlier Detection (Mahalanobis Distance)
library(stats)
library(MASS)
pca <- prcomp(final_data_imputed[, -c(1,2)], center = TRUE, scale. = TRUE) # Exclude ID and Age variables 
pca_data <- as.data.frame(pca$x)
mahalanobis_dist <- mahalanobis(pca_data, colMeans(pca_data), cov(pca_data))
threshold <- qchisq(0.99, df = ncol(pca_data))
outliers <- which(mahalanobis_dist > threshold)

# Exclude outliers
cleaned_data <- final_data_imputed[-outliers, ]
