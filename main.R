# Load methods
source("methods/gauss_seidel.R")
source("methods/Euler-method.R")
source("methods/lagrange_interpolation.R")
source("methods/runge_kutta.R")

# Example: Gauss Seidel
A <- matrix(c(4,-1,0,
              -1,4,-1,
              0,-1,4), nrow=3, byrow=TRUE)

B <- c(15,10,10)

gauss_seidel(A,B)

# Example: Euler Method
f <- function(t,y){
  t + y
}

euler_method(f,10,0,1,1)