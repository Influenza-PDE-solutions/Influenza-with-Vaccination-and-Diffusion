# dss004.R
dss004 <- function(xl, xu, nx, f) {
  fx <- rep(0, nx)
  dx <- (xu - xl)/(nx - 1)
  fx[2:(nx-1)] <- (f[3:nx] - f[1:(nx-2)]) / (2*dx)
  return(fx)
}
