source("R/pca.R")
source("R/eigen.R")
source("R/matrix_ops.R")

# 1. Load data directly
raw <- read.csv("data/cell_metrics.csv", row.names = 1)
neural_signal <- as.matrix(raw[, sapply(raw, is.numeric)])

# 2. Impute missing values with column means
if (any(is.na(neural_signal))) {
  for (j in 1:ncol(neural_signal)) {
    if (any(is.na(neural_signal[, j]))) {
      neural_signal[is.na(neural_signal[, j]), j] <- mean(neural_signal[, j], na.rm = TRUE)
    }
  }
}

# 3. Drop metadata identifier columns
metadata_cols <- c("experiment_container_id", "specimen_id", "tld1_id", "tld2_id", "tlr1_id")
neural_signal <- neural_signal[, !colnames(neural_signal) %in% metadata_cols]

# 4. Remove zero-variance columns (safe now that NAs are imputed)
col_vars <- apply(neural_signal, 2, var)
neural_signal <- neural_signal[, col_vars > 0]

cat(sprintf("Cleaned matrix: %d observations x %d neurons\n",
            nrow(neural_signal), ncol(neural_signal)))

# 5. Scale data (centre and divide by standard deviation)
neural_scaled <- scale(neural_signal)

# 6. Run PCA
pca_result <- compute_pca(neural_scaled, k = 40)
# k = 40 rather than ncol(neural_scaled) = 47
# deflation becomes numerically unstable in the last few near-zero
# components — 40 captures 95% of variance (38 components needed)
# so nothing meaningful is lost

# 7. Print summary metrics
cat(sprintf(
  "\n--- PCA Summary ---\nPC1 explains:               %.1f%%\nPC1+PC2:                    %.1f%%\nk for 80%% variance:         %d components\nk for 95%% variance:         %d components\nTotal components available: %d\n",
  pca_result$var_explained[1] * 100,
  pca_result$cum_var_explained[2] * 100,
  which(pca_result$cum_var_explained >= 0.80)[1],
  which(pca_result$cum_var_explained >= 0.95)[1],
  ncol(neural_scaled)
))


#PC1 top loadings
feat_names <- colnames(neural_scaled)
top_pc1 <- sort(abs(pca_result$loadings[,1]), decreasing=TRUE)[1:5]
names(top_pc1) <- feat_names[order(abs(pca_result$loadings[,1]),
                                   decreasing=TRUE)[1:5]]
print(top_pc1)

# PC2 top loadings
top_pc2 <- order(abs(pca_result$loadings[, 2]), decreasing = TRUE)[1:5]
pc2_named <- pca_result$loadings[top_pc2, 2]
names(pc2_named) <- feat_names[top_pc2]
print(round(pc2_named, 5))
