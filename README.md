# Linear Algebra from Scratch in R

A self-directed study project implementing core linear algebra algorithms from first principles in R. This repository was built as preparation for graduate study in data science to demonstrate a foundational understanding of the mathematics underlying scientific computing.

## What is this project?

Rather than relying on R's highly optimized, built-in matrix operations, this project rebuilds them from the ground up. By coding the mathematical logic manually, this repository bridges the gap between theoretical linear algebra and computational implementation, handling edge cases like floating-point arithmetic and singular matrices.

## What did I build?

I built a custom R package structure with two core mathematical engines, complete with documentation and unit tests:

- **Matrix Operations**: Custom implementations of dot products, matrix multiplication (mat_mul), transposes, and norms—all built without R's %*% or crossprod operators.

- **Determinants & Inverses**: Cofactor expansion for 2x2 and 3x3 matrices, establishing the geometric interpretation of linear transformations.

- **Linear Systems Solver**: A robust Gaussian elimination pipeline (solve_gauss) combining forward elimination (with partial pivoting to prevent zero-division crashes) and back substitution to solve $Ax=b$.

## What did I find?

Building these algorithms exposed the practical challenges of translating exact mathematical definitions into numerically stable implementations.

- **Precision**: `mat_mul()` matches R's `%*%` operator to 10 decimal places across 2×2, 3×3, and non-square matrices.
  
- **Accuracy**: `solve_gauss()` recovers exact solution vectors (e.g., x = (2, 3, −1)) from 3×3 systems with residual ‖Ax − b‖ < 1e-10.
  
- **Stability**: The solver detects singular and near-singular matrices before elimination begins, throwing an informative error rather than returning numerically corrupted results.

- **Reliability**: Validated by a testthat suite of 26 unit tests covering edge cases including zero pivots, mismatched dimensions, singular systems, and identity matrix behaviour — all passing.
  
## Repository Structure

| Folder | Contents |
|--------|----------|
| `R/` | Core implementations: `matrix_ops.R`, `linear_systems.R` |
| `vignettes/` | Mathematical write-ups with derivations and worked examples |
| `tests/` | `testthat` suite validating all algorithmic logic |

## Vignettes

The math behind the code is documented in step-by-step vignettes:

- [01: Vectors, Matrices, and Linear Transformations](vignettes/matrices_01.Rmd)
  
- [02: Solving Linear Systems with Gaussian Elimination](vignettes/linear_systems_02.Rmd)

## Running the Tests

To verify the algorithms locally:

```r
# 1. Install the testing framework
install.packages("testthat")

# 2. Load the core functions into your environment
source("R/matrix_ops.R")
source("R/linear_systems.R")

# 3. Run the specific test files
testthat::test_file("tests/testthat/test_matrix_ops.R")
testthat::test_file("tests/testthat/test_linear_systems.R")
```


