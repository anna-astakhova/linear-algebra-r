#' Tests for eigen.R
#' Run with: source("R/matrix_ops.R"); source("R/eigen.R"); testthat::test_file("tests/testthat/test_eigen.R")

# --- power_iteration ---------------------------

test_that ("power_iterationfinds dominant eigenvalue of 2x2 symmetric matrix", {
  A <- matrix(c(2,1,1,2), nrow= 2)
  result <- power_iteration(A)
  # dominant eigenvalue should be 3
  expect_equal(result$value, 3, tolerance = 1e-10)
})

test_that ("power_iteration eigenvector satisfies Av = λv", {
  A <- matrix(c(4, 2, 2, 3), nrow = 2)
  result <- power_iteration(A)

  v <- result$vector
  lv <- result$value

  # Av - λv should be zero
  residual <- as.vector(A %*% v) - lv*v
  expect_equal(residual, c(0,0), tolerance = 1e-10)
})

test_that("power_iteration returns a unit vector", {
  A <- matrix(c(3, 1, 1, 3), nrow = 2)
  result <- power_iteration(A)

  norm_v <- sqrt(sum(result$vector^2))
  expect_equal(norm_v, 1, tolerance = 1e-10)
})

#--- deflate ---------------------------------------

test_that("deflation removes dominant eigenpair contribution", {
  A <- matrix(c(2, 1, 1, 2), nrow = 2)
  result <- power_iteration(A)
  A2 <- deflate(A, eigenvalue = result$value, v = result$vector)

  # after deflation, the dominant eigenvalue of A2 should be
  # the second eigenvalue of A (which is 1 for this matrix)

  result2 <- power_iteration(A2)
  expect_equal(abs(result2$value), 1, tolerance = 1e-5)
})

#--- eigen_decomp ----------------------------------

test_that("eigen_decomp eigenvalues match R's eigen() for 2x2", {
  A <- matrix(c(4, 2, 2, 3), nrow = 2)

  my  <- eigen_decomp(A)
  ref <- eigen(A)

  # sort both descending before comparing
  expect_equal(sort(my$values, decreasing = TRUE),
               sort(ref$values, decreasing = TRUE),
               tolerance = 1e-5)
})

test_that("eigen_decomp eigenvalues match R's eigen() for 3x3", {
  A <- matrix(c(6, 2, 1, 2, 3, 1, 1, 1, 1), nrow = 3)
  A <- (A + t(A)) / 2    # symmetrise to be safe

  my <- eigen_decomp(A)
  ref <- eigen(A)

  expect_equal(sort(my$values, decreasing = TRUE),
               sort(ref$values, decreasing = TRUE),
               tolerance = 1e-5)
})

test_that("all eigenvectors satisfy Av = λv", {
  A <- matrix(c(5, 2, 2, 5, 1, 0, 1, 0, 3), nrow = 3)
  A <- (A + t(A)) / 2

  result <- eigen_decomp(A)

  for (i in seq_len(ncol(result$vectors))) {
    v <- result$vectors[,i]
    lv <- result$values[i]
    residual <- as.vector( A%*%v) - lv*v
    expect_equal(residual, rep(0,nrow(A)), tolerance = 1e-5)
  }
})

test_that("eigenvectors are orthogonal for symmetric matrix", {
  A <- matrix(c(3, 1, 1, 3), nrow = 2)
  result <- eigen_decomp(A)

  v1 <- result$vectors[,1]
  v2 <- result$vectors[,2]

  # dot product of orthogonal vectors should be 0
  expect_equal(sum(v1*v2), 0, tolerance = 1e-5)
})
