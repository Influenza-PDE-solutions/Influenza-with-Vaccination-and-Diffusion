flu_diffusion <- function(t, u, parms) {
  S = u[1:nx]; V = u[(nx+1):(2*nx)]; E = u[(2*nx+1):(3*nx)]
  I = u[(3*nx+1):(4*nx)]; R = u[(4*nx+1):(5*nx)]
  
  Sx = dss004(xl, xu, nx, S); Sx[c(1, nx)] = 0
  Vx = dss004(xl, xu, nx, V); Vx[c(1, nx)] = 0
  Ex = dss004(xl, xu, nx, E); Ex[c(1, nx)] = 0
  Ix = dss004(xl, xu, nx, I); Ix[c(1, nx)] = 0
  Rx = dss004(xl, xu, nx, R); Rx[c(1, nx)] = 0
  
  Sxx = dss044(xl, xu, nx, S, Sx, 2, 2)
  Vxx = dss044(xl, xu, nx, V, Vx, 2, 2)
  Exx = dss044(xl, xu, nx, E, Ex, 2, 2)
  Ixx = dss044(xl, xu, nx, I, Ix, 2, 2)
  Rxx = dss044(xl, xu, nx, R, Rx, 2, 2)
  
  ut = numeric(5 * nx)
  ut[1:nx]             = d1 * Sxx
  ut[(nx+1):(2*nx)]    = d2 * Vxx
  ut[(2*nx+1):(3*nx)]  = d3 * Exx
  ut[(3*nx+1):(4*nx)]  = d4 * Ixx
  ut[(4*nx+1):(5*nx)]  = d5 * Rxx
  
  return(list(ut))
}
