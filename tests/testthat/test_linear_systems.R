#' Tests for linear_systems.R
#' Run with: ((source("R/linear_systems.R"))); testthat::test_file("tests/testthat/test_linear_systems.R")

# ── row_reduce ───────────────────────────────────────────

test_that("row_reduce produces upper triangular matrux", {
  A <- matrix(c(2,1,-1,-3,-1,2,-2,1,2), nrow=3, byrow=TRUE)
  b <- c(8, -11, -3)
  Ab <- cbind (A,b)
  Ub <- row_reduce(Ab)

  # everything below the diagonal should be (near) zero
  for (i in 2:3) {
    for (j in 1:(i - 1)) {
      expect_equal(as.numeric(Ub[i,j]), 0, tolerance = 1e-10)
    }
  }
})

test_that("row_reduce handles a zero pivot via partial pivoting", {
  # first entry is 0 — would crash without pivoting
  A <- matrix(c(0,2,3,1), nrow=2, byrow=TRUE)
  b <- c(4, 9)
  Ab <- cbind(A, b)
  expect_no_error(row_reduce(Ab))
})

# ── back_substitute ───────────────────────────────────────

test_that("back_substitute solves a simple upper triangular system", {
  Ub <- matrix(c(2,-1,3, 0,1,2), nrow=2, byrow=TRUE)
  x  <- back_substitute(Ub)
  expect_equal(x[2], 2,   tolerance = 1e-10)
  expect_equal(x[1], 2.5, tolerance = 1e-10)
})

# ── solve_gauss ───────────────────────────────────────────

test_that("solve_gauss matches base R solve() on 3x3 system", {
  A <- matrix(c(2,1,-1,-3,-1,2,-2,1,2), nrow=3, byrow=TRUE)
  b <- c(8,-11,-3)
  expect_equal(solve_gauss(A,b), solve(A,b), tolerance = 1e-10)
})

test_that("solve_gauss matches base R solve() on 4x4 system", {
  set.seed(42)
  A <- matrix(sample(1:20, 16), nrow=4)
  b <- sample(1:10, 4)
  expect_equal(solve_gauss(A, b), solve(A, b), tolerance = 1e-8)
})

test_that("solve_gauss solution actually satisfies Ax = b", {
  A <- matrix(c(4,3,6,3), nrow=2, byrow=TRUE)
  b <- c(10, 12)
  x <- solve_gauss(A,b)
  #verify by substituting back: A %*% x should equal b
  expect_equal(as.vector(A %*% x), b, tolerance = 1e-10)
})

test_that("solve_gauss errors on singular matrix", {
  A <- matrix(c(1,2,2,4), nrow=2, byrow=TRUE)  # row 2 = 2 * row 1
  b <- c(3, 6)
  expect_error(solve_gauss(A,b))
})

test_that("solve_gauss errors on non-square matrix", {
  A <- matrix(1:6, nrow=2)
  b <- c(1, 2)
  expect_error(solve_gauss(A, b))
})
