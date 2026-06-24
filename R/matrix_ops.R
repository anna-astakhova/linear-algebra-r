#' Compute the dot product of two vectors
#'
#' @param u Numeric vector
#' @param v Numeric vector of the same length as u
#' @return A single number: the sum of element-wise products

dot_product <- function (u,v) {
  if (length(u) != length(v)) stop ('Vectors must be the same length')
  total <- 0
  for (i in seq_along(u)) {
  total <- total + u[i] * v[i]
  }
  total
}



#' Multiply two matrices from scratch
#'
#' @param A Numeric matrix (m x n)
#' @param B Numeric matrix (n x p)
#' @return Numeric matrix (m x p)


mat_mul <- function (A,B){
  if (ncol(A) != nrow(B)) stop ("Inner dimensions must match: ncol(A) must equal nrow(B)")
    m = nrow (A)
    n = ncol (B)
    o = ncol (A)
    result <- matrix (0, nrow = m, ncol = n)
    for (i in 1:m) {
      for (j in 1:n) {
        for (k in 1:o) {
          result [i, j] <- (A[i,k] * B[k,j]) + result[i,j]
        }
      }
    }
result
}



#' Transpose a matrix from scratch
#'
#' @param A Numeric matrix (m x n)
#' @return Numeric matrix (n x m)
mat_transpose <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  result <- matrix(0, nrow = n, ncol = m)
  for (i in 1:m) {
    for (j in 1:n) {
      result[j, i] <- A[i, j]
    }
  }
  result
}


#' Compute the L2 (Euclidean) norm of a vector
#'
#' @param v Numeric vector
#' @return A single positive number: the length of the vector

vec_norm <- function (V) {
  n <- length (V)
  sum_squares = 0
  for (i in 1:n) {
    sum_squares <- V[i]^2 + sum_squares
  }
  sqrt(sum_squares)
}
