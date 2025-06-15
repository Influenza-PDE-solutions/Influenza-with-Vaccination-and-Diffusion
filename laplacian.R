# Helper: second spatial derivative
laplacian <- function(u, dx) {
  n <- length(u)
  du2 <- numeric(n)
  du2[2:(n-1)] <- (u[3:n] - 2*u[2:(n-1)] + u[1:(n-2)]) / dx^2
  # Neumann BCs (zero flux)
  du2[1] <- du2[2]
  du2[n] <- du2[n-1]
  return(du2)
}