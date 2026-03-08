gauss_seidel <- function(A, B, tol=1e-6){
  n <- nrow(A)
  X <- rep(1,n)

  repeat{

    X_old <- X

    for(i in 1:n){
      cat("Iteration: ", i, " | A[i,]: ", A[i,], " | X: ", X, "\n")
      
      if(i == 1){
        sum1 <- 0
      } else {
        sum1 <- sum(A[i,1:(i-1)] * X[1:(i-1)])
      }
      
      if(i == n){
        sum2 <- 0
      } else {
        sum2 <- sum(A[i,(i+1):n] * X_old[(i+1):n])
      }
      
      X[i] <- (B[i] - sum1 - sum2)/A[i,i]
    }

    if(max(abs(X-X_old)) < tol){
      break
    }
  }

  print(X)
}