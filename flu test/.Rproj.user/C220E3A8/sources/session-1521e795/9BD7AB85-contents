
# Access ODE integrator
library("deSolve");
#
# Access functions for analytical solutions
setwd("E:/R projects/flu test/flu_model");
source("flu_1.R");
source("dss004.R");
source("dss044.R");

#
# Format of output
#
# ip = 1 - graphical (plotted) solutions vs x with
# t as a parameter
#
# ip = 2 - graphical solutions vs t at specific x
#
ip=1;
#
# Grid in x
nx=61;xl=-3;xu=3;
xg=seq(from=xl,to=xu,by=(xu-xl)/(nx-1));
#
# Grid in t
if(ip==1){nout=11;t0=0;tf=60;}
if(ip==2){nout=61;t0=0;tf=60;}
# if(ip==2){nout=61;t0=0;tf=2000;}
tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Parameters
beta=0.5140; betae=0.250; betai=1; betav=0.9;
sigma=1/2; gamma=1/5; delta=1/365;
mu=5.50e-08;
r=1.140e-05; kappa=1.857e-04; alpha=9.30e-06;
theta=1/365;
phi=1/20; d1=0.05; d2=0.05; d3=0.025;
d4=0.001; d5=0;
#
# Display selected parameters
cat(sprintf(
  "\n\n betav = %6.3f phi = %6.3f\n",betav,phi));
#
# ICs
u0=rep(0,5*nx);
for(ix in 1:nx){
  u0[ix] =0.86*exp(-(xg[ix]/1.4)^2);
  u0[ix+nx] =0.10*exp(-(xg[ix]/1.4)^2);
  u0[ix+2*nx]=0;
  u0[ix+3*nx]=0.04*exp(-xg[ix]^2);
  u0[ix+4*nx]=0;
}
ncall=0;
#
# ODE integration
time_taken <- system.time({
  out <- ode(y = u0, times = tout, func = flu_1, parms = NULL)
})
nrow(out)
ncol(out)

#
# Arrays for plotting numerical solutions
S_xplot=matrix(0,nrow=nx,ncol=nout);
V_xplot=matrix(0,nrow=nx,ncol=nout);
E_xplot=matrix(0,nrow=nx,ncol=nout);
I_xplot=matrix(0,nrow=nx,ncol=nout);
R_xplot=matrix(0,nrow=nx,ncol=nout);
for(it in 1:nout){
  for(ix in 1:nx){
    S_xplot[ix,it]=out[it,ix+1];
    V_xplot[ix,it]=out[it,ix+1+nx];
    E_xplot[ix,it]=out[it,ix+1+2*nx];
    I_xplot[ix,it]=out[it,ix+1+3*nx];
    R_xplot[ix,it]=out[it,ix+1+4*nx];
  }
}

#
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
# Calls to ODE routine
cat(sprintf("\n\n ncall = %5d\n\n",ncall));

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
Sm=S_xplot
Vm=V_xplot
Em=E_xplot
Im=I_xplot
Rm=R_xplot
save(Sm, Vm, Em, Im, Rm, file = "MOL_results.RData")

cat("\nTime taken for ODE integration:\n")
print(time_taken)
