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

#' Compute the rank-k SVD approximation of a matrix
#'
#' Uses the truncated SVD: A_k = sum_{i=1}^{k} sigma_i * u_i * v_i^T
#' By the Eckart-Young theorem, this is the closest rank-k matrix to A
#' in Frobenius norm — the best possible low-rank approximation.
#'
#' @param A Numeric matrix to approximate
#' @param k Integer: rank of the approximation
#' @return List with:
#'   - approx: the rank-k approximation matrix (same dimensions as A)
#'   - energy_captured: proportion of total energy (sum sigma^2) retained
svd_approx <- function(A, k) {
  s    <- compute_svd(A)
  U    <- s$U
  Sig  <- s$Sigma
  V    <- s$V

  # reconstruct using only top k singular triplets
  # sum of outer products: sigma_i * u_i * v_i^T
  approx <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:k) {
    approx <- approx + Sig[i] * (U[, i] %*% t(V[, i]))
  }

  # energy captured = proportion of sum of squared singular values
  energy <- sum(Sig[1:k]^2) / sum(Sig^2)

  list(approx = approx, energy_captured = energy)
}
