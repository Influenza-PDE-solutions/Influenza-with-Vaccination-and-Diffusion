flu_reaction <- function(t, u, parms) {
  S = u[1:nx]; V = u[(nx+1):(2*nx)]; E = u[(2*nx+1):(3*nx)]
  I = u[(3*nx+1):(4*nx)]; R = u[(4*nx+1):(5*nx)]
  
  St = Vt = Et = It = Rt = numeric(nx)
  
  for (i in 1:nx) {
    ES = E[i]*S[i]; IS = I[i]*S[i]; EV = E[i]*V[i]; IV = I[i]*V[i]
    IE = I[i]*E[i]; IR = I[i]*R[i]
    
    St[i] = -beta*betae*ES - beta*betai*IS + alpha*IS - phi*S[i] - r*S[i] + delta*R[i] + theta*V[i]
    Vt[i] = -beta*betae*betav*EV - beta*betai*betav*IV + alpha*IV - r*V[i] - theta*V[i] + phi*S[i]
    Et[i] = beta*betae*ES + beta*betai*IS + beta*betae*betav*EV + beta*betai*betav*IV +
      alpha*IE - (r + kappa + sigma)*E[i]
    It[i] = sigma*E[i] - (r + alpha + gamma)*I[i] + alpha*I[i]^2
    Rt[i] = kappa*E[i] + gamma*I[i] - r*R[i] - delta*R[i] + alpha*IR
  }
  
  ut = c(St, Vt, Et, It, Rt)
  return(list(ut))
}

