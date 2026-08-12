#' Compute SVD of a matrix via eigendecomposition of AᵀA
#'
#' Implements A = UΣVᵀ by:
#'   1. Eigendecomposing AᵀA to get V and singular values σ = sqrt(λ)
#'   2. Computing U column by column: uᵢ = Avᵢ / σᵢ
#'
#' This directly exposes the connection between SVD and eigendecomposition.
#' Note: for numerical stability on large matrices, R's built-in svd()
#' uses more sophisticated algorithms — we validate against it below.
#'
#' @param A Numeric matrix (m x n)
#' @return List with: U (m x m), Sigma (vector of singular values), V (n x n)
compute_svd <- function(A) {
  m <- nrow(A)
  n <- ncol(A)

  # step 1: eigendecompose AᵀA — gives right singular vectors V and eigenvalues
  AtA    <- t(A) %*% A
  eig    <- eigen_decomp(AtA, k = min(m, n))

  # singular values: sqrt of eigenvalues (clamp negatives from rounding to 0)
  sigma  <- sqrt(pmax(eig$values, 0))
  V      <- eig$vectors      # n x min(m,n)

  # step 2: compute left singular vectors U = Av / σ
  U <- matrix(0, nrow = m, ncol = length(sigma))
  for (i in seq_along(sigma)) {
    if (sigma[i] > 1e-10) {
      U[, i] <- as.vector(A %*% V[, i]) / sigma[i]
    }
  }

  list(U = U, Sigma = sigma, V = V)
}
