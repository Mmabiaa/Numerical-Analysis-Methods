# Load methods
source("../methods/gauss_seidel.R")
source("../methods/lagrange_interpolation.R")

is_diagonally_dominant <- function(A) {
  n <- nrow(A)
  for (i in 1:n) {
    if (abs(A[i, i]) < sum(abs(A[i, -i]))) {
      return(FALSE)
    }
  }
  return(TRUE)
}

# ---- Example 1: 3x1 - 0.1*x2 - 0.2*x3 = 7.85; 0.1*x1 + 7*x2 - 0.3*x3 = -19.3; 0.3*x1 - 0.2*x2 + 10*x3 = 71.4 ----
cat("\n=============================\n")
cat("Example 1 — Gauss-Seidel Method\n")
cat("System: A*x = B\n")
cat("  3*x1 - 0.1*x2 - 0.2*x3 = 7.85\n")
cat("  0.1*x1 + 7*x2 - 0.3*x3 = -19.3\n")
cat("  0.3*x1 - 0.2*x2 + 10*x3 = 71.4\n")
cat("=============================\n")

A1 <- matrix(c(
  3, -0.1, -0.2,
  0.1, 7, -0.3,
  0.3, -0.2, 10
), nrow = 3, byrow = TRUE)
B1 <- c(7.85, -19.3, 71.4)

cat("Matrix A:\n")
print(A1, digits = 4)
cat("Vector B:", paste(sprintf("%.4f", B1), collapse = " "), "\n\n")

sol1 <- gauss_seidel(A1, B1)
solution1 <- sol1$x
cat("  Solution vector: ", paste(sprintf("%.4f", solution1), collapse = " "), "\n\n")

# ---- Example 4.3: 12*x1 + 3*x2 - 5*x3 = 1; x1 + 5*x2 + 3*x3 = 28; 3*x1 + 7*x2 + 13*x3 = 76 ----
cat("=============================\n")
cat("Example 4.3 — Gauss-Seidel Method\n")
cat("System: A*x = B\n")
cat("  12*x1 + 3*x2 - 5*x3 = 1\n")
cat("  x1 + 5*x2 + 3*x3 = 28\n")
cat("  3*x1 + 7*x2 + 13*x3 = 76\n")
cat("=============================\n")

A2 <- matrix(c(
  12, 3, -5,
  1, 5, 3,
  3, 7, 13
), nrow = 3, byrow = TRUE)
B2 <- c(1, 28, 76)

cat("Matrix A:\n")
print(A2, digits = 4)
cat("Vector B:", paste(sprintf("%.4f", B2), collapse = " "), "\n\n")

sol2 <- gauss_seidel(A2, B2)
solution2 <- sol2$x
cat("  Solution vector: ", paste(sprintf("%.4f", solution2), collapse = " "), "\n\n")

# ---- Example 4.4: Not diagonally dominant (expected to fail) ----
cat("=============================\n")
cat("Example 4.4 — Gauss-Seidel Method\n")
cat("System: A*x = B (matrix is NOT diagonally dominant)\n")
cat("  25*x1 + 5*x2 + x3 = 106.8\n")
cat("  64*x1 + 8*x2 + x3 = 177.2\n")
cat("  144*x1 + 12*x2 + x3 = 279.2\n")
cat("=============================\n")

A3 <- matrix(c(
  25, 5, 1,
  64, 8, 1,
  144, 12, 1
), nrow = 3, byrow = TRUE)
B3 <- c(106.8, 177.2, 279.2)

cat("Matrix A:\n")
print(A3, digits = 4)
cat("Vector B:", paste(sprintf("%.4f", B3), collapse = " "), "\n")
cat("Diagonally dominant? ", is_diagonally_dominant(A3), " (convergence check will fail)\n\n")

tryCatch({
  sol3 <- gauss_seidel(A3, B3)
  solution3 <- sol3$x
  cat("  Solution vector: ", paste(sprintf("%.4f", solution3), collapse = " "), "\n\n")
}, error = function(e) {
  cat("  Result: ", conditionMessage(e), "\n\n")
})

# ---- Example 5.7: Lagrange Interpolating Polynomial (density of motor oil at T = 15 °C) ----
cat("=============================\n")
cat("Example 5.7 — Lagrange Interpolating Polynomial\n")
cat("Evaluate density of unused motor oil at T = 15 °C\n")
cat("=============================\n")

# Data: T (°C) = x, density = f(x)
X_all <- c(0, 20, 40)
Y_all <- c(3.85, 0.800, 0.212)
x_eval <- 15

cat("Data points:\n")
cat("  T (°C)  |  Density\n")
for (i in seq_along(X_all)) {
  cat("  ", sprintf("%.4f", X_all[i]), "   |  ", sprintf("%.4f", Y_all[i]), "\n", sep = "")
}
cat("  Evaluate at T = ", x_eval, " °C\n\n", sep = "")

# First-order (linear): use first two points (0, 3.85) and (20, 0.800)
X1 <- X_all[1:2]
Y1 <- Y_all[1:2]
f1_15 <- lagrange_poly(X1, Y1, x_eval)
cat("First-order Lagrange (linear, 2 points):\n")
cat("  Points: (0, 3.85), (20, 0.800)\n")
cat("  f1(15) = ", sprintf("%.4f", f1_15), "\n\n", sep = "")

# Second-order (quadratic): use all three points
f2_15 <- lagrange_poly(X_all, Y_all, x_eval)
cat("Second-order Lagrange (quadratic, 3 points):\n")
cat("  Points: (0, 3.85), (20, 0.800), (40, 0.212)\n")
cat("  f2(15) = ", sprintf("%.4f", f2_15), "\n\n", sep = "")

cat("Estimated density at T = 15 °C:\n")
cat("  Linear (1st order):  ", sprintf("%.4f", f1_15), "\n", sep = "")
cat("  Quadratic (2nd order): ", sprintf("%.4f", f2_15), "\n", sep = "")
