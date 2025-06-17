# dss044.R
dss044 <- function(xl, xu, nx, f, fx, nl, nu) {
  fxx <- rep(0, nx)
  dx <- (xu - xl)/(nx - 1)
  fxx[2:(nx-1)] <- (f[3:nx] - 2*f[2:(nx-1)] + f[1:(nx-2)]) / (dx^2)
  return(fxx)
}

