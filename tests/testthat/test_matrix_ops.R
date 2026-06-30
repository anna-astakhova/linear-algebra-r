#' Tests for matrix_ops.R
#' Run with: ((source("R/matrix_ops.R"))); testthat::test_file("tests/testthat/test_matrix_ops.R")



# --- dot_product ----------------------------

test_that('dot_product returns correct scalar', {
  expect_equal(dot_product(c(1,2,3), c(4,5,6)), 32)
  expect_equal(dot_product(c(0,1), c(1,0)), 0) # perpendicular -> 0
  expect_equal(dot_product(c(3,4), c(3,4)), 25) # same as norm squared
})

test_that('dot_product errors on mismatched lengths', {
  expect_error(dot_product(c(1,2), c(1,2,3)))
})

#--- mat_mul ---------------------------------

test_that('mat_mul matches base R %*%', {
  A <- matrix (1:4, nrow = 2)
  B <- matrix (5:8, nrow = 2)
  expect_equal(mat_mul(A,B), A %*% B)
})

test_that('mat_mul errors on incompatible dimensions', {
  C <- matrix (1:4, nrow = 2) #2x2
  D <- matrix (1:9, nrow = 3) #3x3
  expect_error(mat_mul(C,D))
})

test_that('mat_mul identity matrix leaves A unchanged', {
  A <- matrix(c(3, 1, 4, 1, 5, 9, 2, 6, 5), nrow = 3)
  I <- diag(3)
  expect_equal(mat_mul(A, I), A)
  expect_equal(mat_mul(I, A), A)
})

#--- mat_transpose ---------------------------

test_that('mat_transpose matches base R t()', {
  A <- matrix (1:6, nrow = 2)
  expect_equal(mat_transpose(A), t(A))
})

test_that('mat_transpose returns original matrix', {
  A <- matrix(sample(1:20), nrow = 4)
  expect_equal(mat_transpose(mat_transpose(A)), A)
})

#--- vec_norm --------------------------------

test_that('vec_norm returns correct Eucledian length', {
  expect_equal(vec_norm(c(3,4)), 5) #classic 3-4-5 triangle
  expect_equal(vec_norm(c(1,0,0)), 1) #unit vector
})

test_that('vec_norm matches base R norm()', {
  v <- c(2.5, 3.1, 7.0)
  expect_equal(vec_norm(v), norm(matrix(v), 'F'))
})

#--- det_2x2 ---------------------------------

test_that('det_2x2 returns correct determinant', {
  expect_equal(det_2x2(matrix(c(3,4,8,6), nrow = 2)), -14)
  expect_equal(det_2x2(matrix(c(2,1,4,2), nrow=2)),  0)  # singular
  expect_equal(det_2x2(diag(2)), 1)                       # identity → det 1
})

test_that('det_2x2 matches base R det()', {
  A <- matrix(c(5, 2, 7, 3), nrow = 2)
  expect_equal(det_2x2(A), det(A))
})

#--- det_3x3 ----------------------------------

test_that('det_3x3 matches base R det()', {
  A <- matrix(c(1,4,7,2,5,8,3,6,10), nrow = 3)
  expect_equal(det_3x3(A), det(A), tolerance = 1e-10)
})

test_that('det_3x3 returns 0 for singular matrix', {
  A <- matrix(1:9, nrow = 3, byrow = TRUE)  # rows are linearly dependent
  expect_equal(det_3x3(A), 0, tolerance = 1e-10)
})

test_that('det_3x3 returns 1 for identity matrix', {
  expect_equal(det_3x3(diag(3)), 1)
})

#--- inv_2x2 ---------------------------------

test_that ('inv_2x2 gives A %*% A_inv = I', {
  A <- matrix(c(3, 4, 8, 6), nrow = 2)
  expect
})

test_that('inv_2x2 mathces base R solve()', {
  A <- matrix(c(4, 3, 3, 2), nrow = 2)
  expect_equal(inv_2x2(A), solve(A), tolerance = 1e-10)
})

test_that('inv_2x2 errors on singular matrix', {
  A <- matrix(c(2, 1, 4, 2), nrow = 2)
  expect_error(inv_2x2(A))
})
