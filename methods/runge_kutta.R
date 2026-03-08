runge_kutta <- function(f,x0,y0,x1,n){

  h <- (x1-x0)/n

  x <- numeric(n+1)
  y <- numeric(n+1)

  x[1] <- x0
  y[1] <- y0

  for(i in 1:n){

    k1 <- h*f(x[i],y[i])
    k2 <- h*f(x[i]+h/2,y[i]+k1/2)
    k3 <- h*f(x[i]+h/2,y[i]+k2/2)
    k4 <- h*f(x[i]+h,y[i]+k3)

    y[i+1] <- y[i] + (k1+2*k2+2*k3+k4)/6
    x[i+1] <- x[i] + h
  }

  plot(x,y,type="l",col="red",main="Runge-Kutta Method")
}