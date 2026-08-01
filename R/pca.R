#' Centre a data matrix by subtracting column means
#'
#' PCA finds directions of maximum variance. Variance is always
#' measured from the mean, so centring is required before computing
#' the covariance matrix. Without it, the first PC would simply
#' point toward the mean of the data, not the direction of spread.
#'
#' @param X Numeric matrix (n x p): n observations, p variables
#' @return Centred matrix of same dimensions: each column has mean 0
centre_data <- function(X) {
  col_means <- colMeans(X)
  # sweep subtracts a vector from each row — equivalent to X - rep(col_means, each=nrow(X))
  sweep(X, 2, col_means, FUN = "-")
}

#' Compute the sample covariance matrix from centred data
#'
#' C = (1/(n-1)) * X_c^T X_c
#' The (i,j) entry of C measures how variables i and j vary together.
#' Diagonal entries are variances; off-diagonal entries are covariances.
#' C is always symmetric and positive semi-definite, guaranteeing
#' real eigenvalues and orthogonal eigenvectors (spectral theorem).
#'
#' @param X_centred Centred numeric matrix (n x p)
#' @return Symmetric covariance matrix (p x p)
covariance_matrix <- function(X_centred) {
  n <- nrow(X_centred)
  (1 / (n - 1)) * t(X_centred) %*% X_centred
}

#' Principal Component Analysis from scratch
#'
#' Implements PCA via eigendecomposition of the sample covariance matrix.
#' This is mathematically equivalent to R's prcomp() but built entirely
#' from our own centre_data(), covariance_matrix(), and eigen_decomp().
#'
#' The k eigenvectors with the largest eigenvalues are the principal
#' components — the orthogonal directions that successively maximise
#' the variance of the projected data.
#'
#' @param X Numeric matrix (n x p): raw (uncentred) data
#' @param k Integer: number of principal components to return (default: all)
#' @return List containing:
#'   - scores:    n x k matrix of PC scores (projected data)
#'   - loadings:  p x k matrix of eigenvectors (principal component directions)
#'   - eigenvalues: numeric vector of length k (variance captured per PC)
#'   - var_explained: proportion of total variance explained by each PC
#'   - cum_var_explained: cumulative proportion of variance explained
compute_pca <- function(X, k = ncol(X)) {
  # step 1: centre
  X_c <- centre_data(X)

  # step 2: covariance matrix
  C <- covariance_matrix(X_c)

  # step 3: eigendecompose — uses your week 5 implementation
  eig <- eigen_decomp(C, k = k)

  # step 4: project data onto top k eigenvectors
  W      <- eig$vectors          # p x k loading matrix
  scores <- X_c %*% W            # n x k score matrix

  # step 5: explained variance
  total_var    <- sum(eigen_decomp(C)$values)
  var_exp      <- eig$values / total_var
  cum_var_exp  <- cumsum(var_exp)

  list(
    scores             = scores,
    loadings           = W,
    eigenvalues        = eig$values,
    var_explained      = var_exp,
    cum_var_explained  = cum_var_exp
  )
}

#' Compare our PCA to R's prcomp() — the ground truth check
#'
#' @param X Raw data matrix
#' @param k Number of components
#' @return Invisible list with comparison results, printed to console
validate_pca <- function(X, k = 2) {
  our   <- compute_pca(X, k)
  ref   <- prcomp(X, center = TRUE, scale. = FALSE)

  cat("=== Eigenvalue comparison ===\n")
  cat("Ours:      ", round(our$eigenvalues, 4), "\n")
  cat("prcomp SD²:", round(ref$sdev[1:k]^2, 4), "\n\n")

  cat("=== Variance explained ===\n")
  cat("Ours:   ", round(our$var_explained[1:k] * 100, 2), "%\n")
  ref_var <- ref$sdev^2 / sum(ref$sdev^2)
  cat("prcomp: ", round(ref_var[1:k] * 100, 2), "%\n\n")

  cat("=== Loadings (may differ in sign — both valid) ===\n")
  cat("Ours PC1:   ", round(our$loadings[,1], 4), "\n")
  cat("prcomp PC1: ", round(ref$rotation[,1], 4), "\n")

  invisible(list(ours = our, ref = ref))
}
