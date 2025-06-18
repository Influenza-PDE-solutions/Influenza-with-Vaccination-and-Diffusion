source("fvm_laplacian.R")
# RHS for FVM system
compute_rhs <- function(S, V, E, I, R, dx) {
  N <- S + V + E + I + R
  ES <- E * S
  IS <- I * S
  EV <- E * V
  IV <- I * V
  IE <- I * E
  IR <- I * R
  
  S_xx <- fvm_laplacian(S, dx)
  V_xx <- fvm_laplacian(V, dx)
  E_xx <- fvm_laplacian(E, dx)
  I_xx <- fvm_laplacian(I, dx)
  R_xx <- fvm_laplacian(R, dx)
  
  dS <- -beta * betae * ES - beta * betai * IS + alpha * IS - phi * S - r * S + delta * R + theta * V + r + d1 * S_xx
  dV <- -beta * betae * betav * EV - beta * betai * betav * IV + alpha * IV - r * V - theta * V + phi * S + d2 * V_xx
  dE <- beta * betae * ES + beta * betai * IS + beta * betae * betav * EV + beta * betai * betav * IV + alpha * IE - (r + kappa + sigma) * E + d3 * E_xx
  dI <- sigma * E - (r + alpha + gamma) * I + alpha * I^2 + d4 * I_xx
  dR <- kappa * E + gamma * I - r * R - delta * R + alpha * IR + d5 * R_xx
  
  list(dS, dV, dE, dI, dR)
}
