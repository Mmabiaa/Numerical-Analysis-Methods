gauss_seidel <- function(A, B, tol = 1e-6, max_iter = 1000, verbose = TRUE) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("Error: Matrix A must be a square matrix")
  }

  n <- nrow(A)
  if (length(B) != n) {
    stop("Error: B must be a vector of length ", n)
  }

  if (!is_diagonally_dominant(A)) {
    stop("Error: Matrix A must be diagonally dominant for convergence")
  }

  X <- rep(1, n)
  k <- 0

  repeat {
    k <- k + 1
    if (k > max_iter) {
      stop("Error: Did not converge within ", max_iter, " iterations")
    }

    X_old <- X

    for (i in 1:n) {
      if (i == 1) {
        sum1 <- 0
      } else {
        sum1 <- sum(A[i, 1:(i - 1)] * X[1:(i - 1)])
      }

      if (i == n) {
        sum2 <- 0
      } else {
        sum2 <- sum(A[i, (i + 1):n] * X_old[(i + 1):n])
      }

      X[i] <- (B[i] - sum1 - sum2) / A[i, i]
    }

    if (max(abs(X - X_old)) < tol) {
      break
    }
  }

  if (verbose) {
    cat("  Converged after ", k, " iteration(s).\n", sep = "")
    cat("  Tolerance used: ", format(tol, scientific = FALSE), "\n", sep = "")
    cat("  Solution:\n")
    for (i in 1:n) {
      cat("    x", i, " = ", sprintf("%.4f", X[i]), "\n", sep = "")
    }
  }

  invisible(list(x = X, iterations = k))
}
