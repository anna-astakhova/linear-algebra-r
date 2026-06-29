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


#' Determinant of a 2x2 matrix
#'
#' Uses the formula: det = ad - bc
#' Geometrically: the factor by which the transformation scales area
#'
#' @param A A 2x2 numeric matrix
#' @return A single number: the determinant
det_2x2 <- function(A) {
  if(!all(dim(A) == c(2,2))) stop ('Matrix must be 2x2')
  A[1,1]*A[2,2] - A[1,2]*A[2,1]
}


#' Determinant of a 3x3 matrix via cofactor expansion
#'
#' Expands along the first column: det = a*M11 - b*M21 + c*M31
#' where Mij is the determinant of the 2x2 submatrix
#' formed by deleting row i and column j
#'
#' @param A A 3x3 numeric matrix
#' @return A single number: the determinant
det_3x3 <- function(A) {
  if(!all(dim(A) == c(3,3))) stop ('Matrix must be 3x3')
  #cofactor along column 1
  A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2]) - A[2,1]*(A[1,2]*A[3,3] - A[1,3]*A[3,2]) + A[3,1] * (A[1,2]*A[2,3]-A[1,3]*A[2,2])
}

#' Inverse of a 2x2 matrix
#'
#' Uses the closed-form formula:
#' A_inv = (1/det(A)) * [[d, -b], [-c, a]]
#'
#' @param A A 2x2 numeric matrix
#' @return The 2x2 inverse matrix
inv_2x2 <- function(A) {
  if(!all(dim(A) == c(2,2))) stop ('Matrix must be 2x2')
  d <- det_2x2(A)
  if(abs(d) < 1e-10) stop ('Matrix is singular (det ≈ 0): no inverse exists')
  (1/d) * matrix(c(A[2,2], -A[1,2], -A[2,1], A[1,1]), nrow = 2)
}

#' Check that A_inv is the true inverse of A
#'
#' Multiplies A by A_inv using mat_mul() and checks
#' the result is close to the identity matrix
#'
#' @param A Original matrix
#' @param A_inv Proposed inverse
#' @param tol Numerical tolerance (default 1e-9)
#' @return TRUE if A %*% A_inv ≈ I, FALSE otherwise
is_inverse <- function(A, A_inv, tol = 1e-9) {
  if (nrow(A) != ncol(A)) stop("Matrix must be square")
  product <- mat_mul(A, A_inv)
  I <- diag(nrow(A))
  all(abs(product - I) < tol)
}
