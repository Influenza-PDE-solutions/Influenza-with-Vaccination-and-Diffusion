# Crank-Nicolson Method for SVEIR Influenza Model in R

ip <- 1 # 1: plot vs x at times t=0,6,...,60 | 2: plot vs t at x=0

# Grid setup
nx <- 61
xl <- -3
xu <- 3
x <- seq(xl, xu, length.out = nx)
dx <- x[2] - x[1]
dx2 <- dx^2

# Time grid
if(ip==1){ nout <- 11; t0 <- 0; tf <- 60 }
if(ip==2){ nout <- 61; t0 <- 0; tf <- 60 }
tout <- seq(from=t0, to=tf, length.out=nout)

dt <- 0.01
nt <- length(seq(t0, tf, by=dt))
time <- seq(t0, tf, length.out = nt)

# Parameters
beta <- 0.5140
betae <- 0.250
betai <- 1
betav <- 0.9
sigma <- 1/2
gamma <- 1/5
delta <- 1/365
r <- 1.140e-05
kappa <- 1.857e-04
alpha <- 9.30e-06
theta <- 1/365
phi <- 1/20
mu <- 5.50e-08
d1 <- 0.05
d2 <- 0.05
d3 <- 0.025
d4 <- 0.001
d5 <- 0.0

# Initial conditions
S <- matrix(0, nx, nt)
V <- matrix(0, nx, nt)
E <- matrix(0, nx, nt)
I <- matrix(0, nx, nt)
R <- matrix(0, nx, nt)

S[,1] <- 0.86 * exp(-(x/1.4)^2)
V[,1] <- 0.10 * exp(-(x/1.4)^2)
I[,1] <- 0.04 * exp(-x^2)

# Crank-Nicolson matrix
create_cn_matrix <- function(d, nx, dt, dx2){
  lambda <- d * dt / (2 * dx2)
  main <- rep(1 + 2 * lambda, nx)
  off <- rep(-lambda, nx - 1)
  A <- diag(main)
  for(i in 1:(nx-1)){
    A[i,i+1] <- off[i]
    A[i+1,i] <- off[i]
  }
  return(A)
}

A_S <- create_cn_matrix(d1, nx, dt, dx2)
A_V <- create_cn_matrix(d2, nx, dt, dx2)
A_E <- create_cn_matrix(d3, nx, dt, dx2)
A_I <- create_cn_matrix(d4, nx, dt, dx2)
A_R <- create_cn_matrix(d5, nx, dt, dx2)

solve_cn <- function(A, rhs){
  return(solve(A, rhs))
}

# Time stepping
for(k in 1:(nt-1)){
  St <- S[,k]; Vt <- V[,k]; Et <- E[,k]; It <- I[,k]; Rt <- R[,k]
  ES <- Et * St; IS <- It * St; EV <- Et * Vt; IV <- It * Vt
  IE <- It * Et; IR <- It * Rt
  
  F_S <- -beta*betae*ES - beta*betai*IS + alpha*IS - phi*St - r*St + delta*Rt + theta*Vt
  F_V <- -beta*betae*betav*EV - beta*betai*betav*IV + alpha*IV - r*Vt - theta*Vt + phi*St
  F_E <- beta*betae*ES + beta*betai*IS + beta*betae*betav*EV + beta*betai*betav*IV +
    alpha*IE - (r + kappa + sigma)*Et
  F_I <- sigma*Et - (r + alpha + gamma)*It + alpha * It^2
  F_R <- kappa*Et + gamma*It - r*Rt - delta*Rt + alpha*IR
  
  S[,k+1] <- solve_cn(A_S, St + dt*F_S)
  V[,k+1] <- solve_cn(A_V, Vt + dt*F_V)
  E[,k+1] <- solve_cn(A_E, Et + dt*F_E)
  I[,k+1] <- solve_cn(A_I, It + dt*F_I)
  R[,k+1] <- solve_cn(A_R, Rt + dt*F_R)
}

# Prepare output matrices
S_xplot <- matrix(0, nx, nout)
V_xplot <- matrix(0, nx, nout)
E_xplot <- matrix(0, nx, nout)
I_xplot <- matrix(0, nx, nout)
R_xplot <- matrix(0, nx, nout)

time_indices <- round(seq(1, nt, length.out=nout))

for(j in 1:nout){
  S_xplot[,j] <- S[,time_indices[j]]
  V_xplot[,j] <- V[,time_indices[j]]
  E_xplot[,j] <- E[,time_indices[j]]
  I_xplot[,j] <- I[,time_indices[j]]
  R_xplot[,j] <- R[,time_indices[j]]
}
# Display numerical solutions (for t = 0, 60)
if(ip==1){
  for(it in 1:nout){
    if((it-1)*(it-11)==0){
      cat(sprintf("\n\n t x S(x,t)
V(x,t)"));
      cat(sprintf("\n E(x,t) I(x,t)
R(x,t)"));
      for(ix in 1:nx){
        cat(sprintf("\n %6.1f%7.2f%12.5f%12.5f",
                    tout[it],xg[ix],S_xplot[ix,it],V_xplot[ix,it]));
        cat(sprintf("\n%14.5f%12.5f%12.5f",
                    E_xplot[ix,it],I_xplot[ix,it],R_xplot[ix,it]));
      }
    }
  }
}
if(ip==2){
  for(it in 1:nout){
    if((it-1)*(it-61)==0){
      cat(sprintf("\n\n t x S(x,t)
V(x,t)"));
      cat(sprintf("\n E(x,t) I(x,t)
R(x,t)"));
      for(ix in 1:nx){
        cat(sprintf("\n %6.1f%7.2f%12.5f%12.5f",
                    tout[it],xg[ix],S_xplot[ix,it],V_xplot[ix,it]));
        cat(sprintf("\n%14.5f%12.5f%12.5f",
                    E_xplot[ix,it],I_xplot[ix,it],R_xplot[ix,it]));
      }
    }
  }
}
#






# Plot S,V,E,I,R numerical solutions
#
# vs x with t as a parameter, t = 0,6,...,60
if(ip==1){
  par(mfrow=c(1,1));
  matplot(x=xg,y=S_xplot,type="l",xlab="x",
          ylab="S(x,t), t=0,6,...,60",xlim=c(xl,xu),
          lty=1,main="S(x,t); t=0,6,...,60;",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=xg,y=V_xplot,type="l",xlab="x",
          ylab="V(x,t), t=0,6,...,60",xlim=c(xl,xu),
          lty=1,main="V(x,t); t=0,6,...,60;",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=xg,y=E_xplot,type="l",xlab="x",
          ylab="E(x,t), t=0,6,...,60",xlim=c(xl,xu),
          lty=1,main="E(x,t); t=0,6,...,60;",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=xg,y=I_xplot,type="l",xlab="x",
          ylab="I(x,t), t=0,6,...,60",xlim=c(xl,xu),
          lty=1,main="I(x,t); t=0,6,...,60;",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=xg,y=R_xplot,type="l",xlab="x",
          ylab="R(x,t), t=0,6,...,60",xlim=c(xl,xu),
          lty=1,main="R(x,t); t=0,6,...,60;",lwd=2);
}
#
# vs t at x = 0, t = 0,1,...,60
if(ip==2){
  S_tplot=rep(0,nout);V_tplot=rep(0,nout);
  E_tplot=rep(0,nout);I_tplot=rep(0,nout);
  R_tplot=rep(0,nout);
  for(it in 1:nout){
    S_tplot[it]=S_xplot[31,it];
    V_tplot[it]=V_xplot[31,it];
    E_tplot[it]=E_xplot[31,it];
    I_tplot[it]=I_xplot[31,it];
    R_tplot[it]=R_xplot[31,it];
  }
  par(mfrow=c(1,1));
  matplot(x=tout,y=S_tplot,type="l",xlab="t",
          ylab="S(x,t), x = 0",xlim=c(t0,tf),lty=1,
          main="S(x,t); x = 0",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=tout,y=V_tplot,type="l",xlab="t",
          ylab="V(x,t), x = 0",xlim=c(t0,tf),lty=1,
          main="V(x,t); x = 0",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=tout,y=E_tplot,type="l",xlab="t",
          ylab="E(x,t), x = 0",xlim=c(t0,tf),lty=1,
          main="E(x,t); x = 0",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=tout,y=I_tplot,type="l",xlab="t",
          ylab="I(x,t), x = 0",xlim=c(t0,tf),lty=1,
          main="I(x,t); x = 0",lwd=2);
  par(mfrow=c(1,1));
  matplot(x=tout,y=R_tplot,type="l",xlab="t",
          ylab="R(x,t), x = 0",xlim=c(t0,tf),lty=1,
          main="R(x,t); x = 0",lwd=2);
}
Sc=S_xplot
Vc=V_xplot
Ec=E_xplot
Ic=I_xplot
Rc=R_xplot
save(Sc, Vc, Ec, Ic, Rc, file = "Crank-Nicolson_results.RData")



