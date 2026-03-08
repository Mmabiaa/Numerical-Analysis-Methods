euler_method <- function(f,n,t0,t1,y0){

  h <- (t1-t0)/n

  t <- numeric(n+1)
  y <- numeric(n+1)

  t[1] <- t0
  y[1] <- y0

  for(i in 1:n){
    t[i+1] <- t[i] + h
    y[i+1] <- y[i] + h*f(t[i],y[i])
  }

  plot(t,y,type="l",main="Euler Method")
}