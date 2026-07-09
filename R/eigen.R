#' Power iteration to find the dominant eigenvector
#'
#' Repeatedly applies matrix A to a random vector and normalizes it.
#' The vector converges to the eigenvector associated with the largest
#' absolute eigenvalue. The eigenvalue is then recovered via the
#' Rayleigh quotient: λ = v^T A v.
#'
#' @param A A square numeric matrix.
#' @param n_iter Maximum number of iterations to run (default 1000).
#' @param tol Convergence tolerance; loop breaks when the vector changes by less than this amount (default 1e-10).
#' @return A list containing: \code{vector} (the dominant eigenvector) and \code{value} (the dominant eigenvalue).
#' @export
power_iteration <- function(A, n_iter = 1000, tol = 1e-10) {
  n <- nrow(A)
  v <- rnorm(n)
  v <- v / vec_norm(v)

  for (i in seq_len(n_iter)) {
    v_new <- A %*% v
    v_new <- as.vector(v_new)
    v_new <- v_new / vec_norm(v_new)

    if (vec_norm(v_new - v) < tol) {
      v <- v_new
      break
    }
    v <- v_new
  }

  eigenvalue <- as.numeric(t(v) %*% A %*% v)
  list(vector = v, value = eigenvalue)
}
