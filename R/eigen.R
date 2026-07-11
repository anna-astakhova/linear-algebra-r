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
  # start with a random unit vector
  v <- rnorm(n)
  v <- v / vec_norm(v) # normalise using vec_norm custom function


  for (i in seq_len(n_iter)) {
    v_new <- A %*% v # apply the transformation
    v_new <- as.vector(v_new)
    v_new <- v_new / vec_norm(v_new)  # normalise

    # check convergence: did the vector stop changing?
    if (vec_norm(v_new - v) < tol) {
      v <- v_new
      break
    }
    v <- v_new
  }

  # recover eigenvalue via Rayleigh quotient: λ = v'Av
  eigenvalue <- as.numeric(t(v) %*% A %*% v)
  list(vector = v, value = eigenvalue)
}

#' Deflate a matrix by removing a known eigenpair
#'
#' Computes A_deflated = A - lambda * (v %*% t(v)).
#' After deflation, running power iteration on A_deflated
#' finds the next largest eigenvector.
#'
#' @param A Square numeric matrix.
#' @param eigenvalue Scalar: the known dominant eigenvalue.
#' @param v Numeric vector: the known dominant eigenvector.
#' @return Deflated matrix of the same dimensions as A.
#' @export
deflate <- function(A, eigenvalue, v) {
  # Calculate the outer product and subtract the weighted footprint
  A - eigenvalue * (v%*%t(v))
}

#' Eigendecomposition via power iteration and deflation
#'
#' Finds the top k eigenpairs of a symmetric matrix by repeatedly
#' applying power iteration and deflating the result.
#'
#' @param A Square symmetric numeric matrix.
#' @param k Number of eigenpairs to find (default: all, i.e., nrow(A)).
#' @return List with: \code{values} (numeric vector of k eigenvalues, descending),
#'         and \code{vectors} (n x k matrix of corresponding eigenvectors as columns).
#' @export
eigen_decomp <- function(A, k = nrow(A)) {
  n <- nrow(A)

  # Pre-allocate containers for maximum speed
  eigenvalues  <- numeric(k)
  eigenvectors <- matrix(0, nrow = n, ncol = k)

  # Create a working copy so we don't destroy the original matrix
  A_current <- A

  for (i in seq_len(k)) {
    # Find the dominant current in the working matrix
    result <- power_iteration(A_current)

    # Save the value and vector into our pre-allocated containers
    eigenvalues[i]    <- result$value
    eigenvectors[, i] <- result$vector

    # Mathematically remove this current from the working matrix
    A_current <- deflate(A_current, result$value, result$vector)
  }

  list(values = eigenvalues, vectors = eigenvectors)
}
