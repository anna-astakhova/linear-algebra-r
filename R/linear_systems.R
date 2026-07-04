#' Forward elimination with partial pivoting
#'
#' Reduces the augmented matrix [A | b] to upper triangular form [U | b'].
#' Mathematically, this applies elementary row operations. For each pivot column k
#' (from 1 to n-1), it finds the row i >= k where abs(A[i,k]) is maximized and swaps
#' row i with row k (partial pivoting to prevent division by zero).
#' Then, for all rows j > k, it updates row j by subtracting a scaled version of row k:
#' Row_j = Row_j - (A[j,k] / A[k,k]) * Row_k
#' This systematically forces all entries below the main diagonal to become zero.
#'
#' @param Ab Augmented matrix [A | b], dimension n x (n+1)
#' @return Upper triangular augmented matrix [U | b']
row_reduce <- function(Ab) {
  n <- nrow(Ab)

  for (col in 1:(n-1)) {

    # --- partial pivoting: find row with largest absolute value in pivot col
    max_row <- which.max(abs(Ab[col:n, col])) + col - 1

    # swap current row with max_row if they differ
    if (max_row != col) {
      Ab[c(max_row, col), ] <- Ab[c(col, max_row), ]
    }

    # guard: if pivot is still 0 after swapping, matrix is singular
    if (abs(Ab[col,col]) < 1e-10) {
      stop ("Matrix is singular or nearly singular: system may have no unique solution")
    }

    # --- eliminate entries below the pivot
    for (row in (col+1):n) {
      multiplier <- Ab[row,col]/Ab[col,col]
      Ab[row, ] <- Ab[row, ] - multiplier * Ab[col, ]
    }
  }
  # --- FINAL GUARD: check the very last diagonal element
  if (abs(Ab[n, n]) < 1e-10) {
    stop("Matrix is singular or nearly singular: system may have no unique solution")
  }
  Ab
}

#' Back substitution on an upper triangular augmented matrix
#'
#' Mathematically solves the upper triangular system Ux = b' from the bottom up.
#' For the final variable, it computes x_n = b'_n / U_nn.
#' For all preceding variables i (from n-1 down to 1), it isolates x_i using
#' the already-computed variables:
#' x_i = (b'_i - sum_{j=i+1}^n (U_ij * x_j)) / U_ii
#'
#' @param Ub Upper triangular augmented matrix [U | b'], dimension n x (n+1)
#' @return Numeric vector x of length n: the solution to the system
back_substitute <- function(Ub) {
  n <- nrow(Ub)
  x <- numeric(n)  # initialise solution vector with zeros

  # --- loop backwards from the bottom row up to the top
  for (i in n:1) {

    # isolate the constant on the right-hand side and the diagonal coefficient
    rhs <- Ub[i, n + 1]
    pivot <- Ub[i, i]

    # --- calculate contributions from already-solved variables
    # guard: the bottom row has no previously solved variables to substitute
    if (i == n) {
      substitution <- 0
    } else {
      # sum the products of known coefficients and previously solved x values
      substitution <- sum(Ub[i, (i + 1):n] * x[(i + 1):n])
    }

    # --- solve for the current variable
    x[i] <- (rhs - substitution) / pivot
  }

  return(x)
}

#' Solve a linear system Ax = b using Gaussian elimination
#'
#' Mathematically solves for the vector x in the matrix equation Ax = b.
#' It combines forward elimination (transforming A into an upper triangular
#' matrix U while applying the same operations to b) and back substitution
#' (solving the resulting Ux = b' iteratively).
#'
#' @param A Square numeric matrix (n x n): coefficients
#' @param b Numeric vector of length n: right-hand side
#' @return Numeric vector x of length n satisfying Ax = b
solve_gauss <- function(A, b) {
  n <- nrow(A)
  if (nrow(A) != ncol(A)) stop("A must be a square matrix")
  if (length(b) != n)     stop("b must have the same length as nrow(A)")

  # build augmented matrix [A | b]
  Ab <- cbind(A, b)

  # forward elimination
  Ub <- row_reduce(Ab)

  # back substitution
  back_substitute(Ub)
}
