#finit volume method
source("compute_rhs.R")

# ip = 1 - graphical (plotted) solutions vs x with
# t as a parameter
#
# ip = 2 - graphical solutions vs t at specific x
#
ip=1;
# Grid in x
nx=61;xl=-3;xu=3;
xg=seq(from=xl,to=xu,by=(xu-xl)/(nx-1));
#
# Grid in t
if(ip==1){nout=11;t0=0;tf=60;}
if(ip==2){nout=61;t0=0;tf=60;}
# if(ip==2){nout=61;t0=0;tf=2000;}
tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));


# PARAMETERS
beta  <- 0.5140
betae <- 0.250
betai <- 1
betav <- 0.9
sigma <- 1/2
gamma <- 1/5
delta <- 1/365
r     <- 1.140e-05
kappa <- 1.857e-04
alpha <- 9.30e-06
theta <- 1/365
phi   <- 1/20
mu    <- 5.50e-08
d1 <- 0.05
d2 <- 0.05
d3 <- 0.025
d4 <- 0.001
d5 <- 0.0

# SPACE-TIME GRID

x <- seq(xl, xu, length.out = nx)
dx <- x[2] - x[1]

dt <- 0.01
nt <- floor((tf - t0)/dt) + 1
time <- seq(t0, tf, length.out = nt)

# INITIAL CONDITIONS
S <- matrix(0, nx, nt)
V <- matrix(0, nx, nt)
E <- matrix(0, nx, nt)
I <- matrix(0, nx, nt)
R <- matrix(0, nx, nt)

S[, 1] <- 0.86 * exp(-(x / 1.4)^2)
V[, 1] <- 0.10 * exp(-(x / 1.4)^2)
I[, 1] <- 0.04 * exp(-x^2)

time_taken <- system.time({
# RK4 TIME INTEGRATION
for (t in 1:(nt - 1)) {
  S_t <- S[, t]; V_t <- V[, t]; E_t <- E[, t]; I_t <- I[, t]; R_t <- R[, t]
  
  k1 <- compute_rhs(S_t, V_t, E_t, I_t, R_t, dx)
  k2 <- compute_rhs(S_t + 0.5 * dt * k1[[1]],
                    V_t + 0.5 * dt * k1[[2]],
                    E_t + 0.5 * dt * k1[[3]],
                    I_t + 0.5 * dt * k1[[4]],
                    R_t + 0.5 * dt * k1[[5]], dx)
  k3 <- compute_rhs(S_t + 0.5 * dt * k2[[1]],
                    V_t + 0.5 * dt * k2[[2]],
                    E_t + 0.5 * dt * k2[[3]],
                    I_t + 0.5 * dt * k2[[4]],
                    R_t + 0.5 * dt * k2[[5]], dx)
  k4 <- compute_rhs(S_t + dt * k3[[1]],
                    V_t + dt * k3[[2]],
                    E_t + dt * k3[[3]],
                    I_t + dt * k3[[4]],
                    R_t + dt * k3[[5]], dx)
  
  S[, t+1] <- S_t + (dt/6) * (k1[[1]] + 2*k2[[1]] + 2*k3[[1]] + k4[[1]])
  V[, t+1] <- V_t + (dt/6) * (k1[[2]] + 2*k2[[2]] + 2*k3[[2]] + k4[[2]])
  E[, t+1] <- E_t + (dt/6) * (k1[[3]] + 2*k2[[3]] + 2*k3[[3]] + k4[[3]])
  I[, t+1] <- I_t + (dt/6) * (k1[[4]] + 2*k2[[4]] + 2*k3[[4]] + k4[[4]])
  R[, t+1] <- R_t + (dt/6) * (k1[[5]] + 2*k2[[5]] + 2*k3[[5]] + k4[[5]])
}
})
# Find time indices corresponding to tout
it_index <- sapply(tout, function(tval) which.min(abs(time - tval)))

# Extract solution slices for display
S_xplot <- S[, it_index]
V_xplot <- V[, it_index]
E_xplot <- E[, it_index]
I_xplot <- I[, it_index]
R_xplot <- R[, it_index]


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

Sf=S_xplot
Vf=V_xplot
Ef=E_xplot
If=I_xplot
Rf=R_xplot
save(Sf, Vf, Ef, If, Rf, file = "FVM_results.RData")

cat("\nTime taken for ODE integration:\n")
print(time_taken)