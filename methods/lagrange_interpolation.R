lagrange_poly <- function(X,Y,x){

  n <- length(X)
  result <- 0

  for(i in 1:n){

    term <- Y[i]

    for(j in 1:n){

      if(j != i){
        term <- term * (x - X[j])/(X[i] - X[j])
      }

    }

    result <- result + term

  }

  return(result)
}