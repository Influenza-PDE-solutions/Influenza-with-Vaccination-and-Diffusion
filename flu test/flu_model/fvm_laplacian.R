# FVM Laplacian: Robust and vectorized
fvm_laplacian <- function(u, dx) {
  n <- length(u)
  flux <- numeric(n)
  
  # Interior points
  flux[2:(n - 1)] <- (u[3:n] - 2 * u[2:(n - 1)] + u[1:(n - 2)]) / dx^2
  
  # Neumann boundary (zero flux): use mirrored values
  flux[1] <- (u[2] - u[1]) / dx^2  # ∂u/∂x = 0 at left
  flux[n] <- (u[n - 1] - u[n]) / dx^2  # ∂u/∂x = 0 at right
  
  return(flux)
}