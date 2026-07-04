#' Compute the dot product of two vectors
#'
#' Mathematically calculates the scalar product of two vectors u and v.
#' It computes the sum of the products of their corresponding entries:
#' u . v = sum_{i=1}^n (u_i * v_i)
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
#' Mathematically computes the matrix product C = AB.
#' If A is an (m x o) matrix and B is an (o x n) matrix, the resulting
#' matrix C has dimensions (m x n). Each element C_ij is the dot product
#' of the i-th row of A and the j-th column of B:
#' C_ij = sum_{k=1}^o (A_ik * B_kj)
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
#' Mathematically reflects a matrix over its main diagonal.
#' For an original matrix A of dimensions (m x n), the transpose A^T
#' has dimensions (n x m) where the element at the i-th row and j-th
#' column of A^T equals the element at the j-th row and i-th column of A:
#' (A^T)_ij = A_ji
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
#' Mathematically calculates the geometric length (magnitude) of a vector v.
#' It uses the Pythagorean theorem extended to n dimensions:
#' ||v||_2 = sqrt(sum_{i=1}^n (v_i^2))
#'
#' @param V Numeric vector
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
#' Mathematically computes the determinant for a 2x2 matrix A = [[a, b], [c, d]].
#' det(A) = (a * d) - (b * c).
#' Geometrically, this represents the oriented area scaling factor of the
#' linear transformation described by the matrix.
#'
#' @param A A 2x2 numeric matrix
#' @return A single number: the determinant
det_2x2 <- function(A) {
  if(!all(dim(A) == c(2,2))) stop ('Matrix must be 2x2')
  A[1,1]*A[2,2] - A[1,2]*A[2,1]
}

#' Determinant of a 3x3 matrix via cofactor expansion
#'
#' Mathematically computes the determinant using Laplace expansion down the
#' first column. For a 3x3 matrix, it breaks the calculation into three 2x2
#' determinants (minors) multiplied by alternating signs:
#' det(A) = A_11*(A_22*A_33 - A_23*A_32)
#'        - A_21*(A_12*A_33 - A_13*A_32)
#'        + A_31*(A_12*A_23 - A_13*A_22)
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
#' Mathematically calculates the matrix A^(-1) such that A * A^(-1) = I.
#' For a 2x2 matrix [[a, b], [c, d]], it uses the closed-form formula:
#' A^(-1) = (1 / det(A)) * [[d, -b],
#'                          [-c, a]]
#' If the determinant is 0, the matrix is singular and cannot be inverted.
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
#' Mathematically verifies the identity A * A^(-1) = I.
#' It multiplies the original matrix A by the proposed inverse A_inv and
#' checks if the resulting (n x n) matrix equals the Identity matrix I
#' (1s on the diagonal, 0s elsewhere) within a small numerical tolerance
#' to account for floating-point arithmetic errors.
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
