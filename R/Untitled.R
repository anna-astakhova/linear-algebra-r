usethis::use_testthat()
usethis::use_vignette('matrices_01')
install.packages(c("usethis", "devtools", "testthat", "knitr"))
install.packages(c("rmarkdown", "ggplot2", "imager"))

usethis::use_r("matrices")

m <- matrix(1:6, nrow = 2, byrow = TRUE)
for (i in 1:2) {
  for (j in 1:3) {
    print(m[i,j])
  }
}


a <- matrix (1:4, nrow = 2, byrow = TRUE)
b <- matrix (5:8, nrow = 2, byrow = TRUE)
c <- array (c(a,b), dim = c(2,2,2))
print(c)


for (l in 1:2) {
  for (i in 1:2) {
    for (j in 1:2) {
      print(c[i,j,l])
    }
  }
}


print(a)

for (columns in 1:ncol(a)) {
  print(a[columns])
}




f <- matrix (1:4, nrow = 2, byrow = TRUE)
u <- matrix (5:8, nrow = 2, byrow = TRUE)
a <- matrix (0, nrow = 2, ncol = 2)
for (i in 1:2) {
  for (j in 1:2) {
   a[i,j] = f[i,j] +u[i,j]
  }
}
print(a)

A<- matrix (1:4, nrow = 2, byrow = TRUE)
B<- matrix (5:8, nrow = 2, byrow = TRUE)

C <- matrix (0, nrow = 2, byrow = TRUE)
matmul <- function (A,B) {
if (ncol(A) != nrow(B)) stop ('wrong dimensions')
  m <- nrow(A)
  n <- ncol(B)
  o <- ncol(A)
  result <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      for (k in 1:o) {
        result[i,j] = (A[i,k] * B [k,j]) + result[i,j]
      }
    }
  }
return(result)
}

# 1. Define your inputs
A <- matrix(1:4, nrow = 2, byrow = TRUE)
B <- matrix(5:8, nrow = 2, byrow = TRUE)

# 2. Run your beautiful function and capture the output
C <- matmul(A, B)

# 3. Print the result
print(C)



transpose <- function (A) {
  if (nrow(A) != ncol(A)) stop ('wrong dimensions')
  result <- matrix(0, nrow = n, ncol = m)
  m <- ncols(A)
  n <- nrows (A)
  for (i in 1:m) {
    for (j in 1:n) {
      result[j,i] = A[i,j]
    }
  }
}


vec_norm <- function (V) {
  n <- length (V)
  sum_squares = 0
  for (i in 1:n) {
      sum_squares = V[i]^2 + sum_squares
  }
  return(sqrt(sum_squares))
}

a<- 1:10
b<- vec_norm(a)

my_vector <- c(3, 4)
vec_norm(my_vector)

det_2x2 <- function(A) {
  if(!all(dim(A) == c(2,2))) stop ('not 2x2 matrix')
  A[1,1]*A[2,2] - A[1,2]*A[2,1]
}

A <- matrix (1:4, nrow =2 , byrow = TRUE)
det_2x2(A)

inverse <- function(A) {
  if(!all(dim(A) == c(2,2))) stop ('mm')
  d <- det_2x2(A)
  if(abs(d) < 1e-10) stop ('nx')
  (1/d) * matrix(c(A[2,2], -A[1,2], -A[2,1], A[1,1]))
}
inverse(A)

is_inverse <- function(A, A_inv, tol = 1e-9) {
  if (nrow(A) != ncol(A)) stop("Matrix must be square")
  product <- mat_mul(A, A_inv)
  I <- diag(nrow(A))
  all(abs(product - I) < tol)
}


row_reduce <- function (Ab) {
  n <- nrow(Ab)
  for (col in 1:(n-1)) {
  max_row <- which.max(abs(Ab[col:n,col])) + col - 1
  if (max_row != col) { Ab[c(max_row, col),] <- Ab[c(col, max_row),]
   }
if (abs(Ab[col,col]) < 1e-10 ) {
  stop ('Matrix is singular') }
for (row in (col+1):n) {
  multiplier <- Ab[row,col]/Ab[col,col]
  Ab[row,] <- Ab[row,] - (multiplier * Ab[col,])
}
  }
  return(Ab)
}

back_substitution <- function (Ub) {
  n<- nrow(Ub)
  x<- numeric(n)
  for (i in n:1) {
    rhs <- Ub [i, n+1]
    pivot <- Ub [i,i]
    substitution <- sum(Ub[i, (i+1):n] * x[(i+1):n])
    x[i]=(rhs-substitution)/pivot
  }
  return(x)
}

file.create(file.path("tests", "testthat", "test_linear_systems.R"))
