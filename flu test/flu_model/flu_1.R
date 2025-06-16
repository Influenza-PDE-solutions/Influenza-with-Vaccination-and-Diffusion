flu_1=function(t,u,parms){
  #
  # Function flu_1 computes the t derivative vector
  # of the S,V,E,I,R vectors
  #
  # One vector to five vectors
  S=rep(0,nx);V=rep(0,nx);
  E=rep(0,nx);I=rep(0,nx);
  R=rep(0,nx);
  for(i in 1:nx){
    S[i]=u[i];
    V[i]=u[i+nx];
    E[i]=u[i+2*nx];
    I[i]=u[i+3*nx];
    R[i]=u[i+4*nx];
  }
  #
  # Boundary conditions
  Sx=dss004(xl,xu,nx,S);
  Vx=dss004(xl,xu,nx,V);
  Ex=dss004(xl,xu,nx,E);
  Ix=dss004(xl,xu,nx,I);
  Rx=dss004(xl,xu,nx,R);
  Sx[1]=0;Sx[nx]=0;
  Vx[1]=0;Vx[nx]=0;
  Ex[1]=0;Ex[nx]=0;
  Ix[1]=0;Ix[nx]=0;
  Rx[1]=0;Rx[nx]=0;
  nl=2;nu=2;
  #
  # Sxx to Rxx
  Sxx=dss044(xl,xu,nx,S,Sx,nl,nu);
  Vxx=dss044(xl,xu,nx,V,Vx,nl,nu);
  Exx=dss044(xl,xu,nx,E,Ex,nl,nu);
  Ixx=dss044(xl,xu,nx,I,Ix,nl,nu);
  Rxx=dss044(xl,xu,nx,R,Rx,nl,nu);
  #
  # PDEs
  b=beta;be=betae;bi=betai;bv=betav;
  a=alpha;p=phi;d=delta;t=theta;k=kappa;
  s=sigma;g=gamma;
  St=rep(0,nx);Vt=rep(0,nx);
  Et=rep(0,nx);It=rep(0,nx);
  Rt=rep(0,nx);
  for(i in 1:nx){
    ES=E[i]*S[i];
    IS=I[i]*S[i];
    EV=E[i]*V[i];
    IV=I[i]*V[i];
    IE=I[i]*E[i];
    IR=I[i]*R[i];
    St[i]=-b*be*ES-b*bi*IS+a*IS-p*S[i]-r*S[i]
    +d*R[i]+t*V[i]+r+d1*Sxx[i];
    Vt[i]=-b*be*bv*EV-b*bi*bv*IV+a*IV-r*V[i]
    -t*V[i]+p*S[i]+d2*Vxx[i];
    Et[i]=b*be*ES+b*bi*IS+b*be*bv*EV+b*bi*bv*IV
    +a*IE-(r+k+s)*E[i]+d3*Exx[i];
    It[i]=s*E[i]-(r+a+g)*I[i]+a*I[i]^2+d4*Ixx[i];
    Rt[i]=k*E[i]+g*I[i]-r*R[i]-d*R[i]+a*IR+d5*Rxx[i];
  }
  #
  # Five vectors to one vector
  ut=rep(0,5*nx);
  for(i in 1:nx){
    ut[i] =St[i];
    ut[i+nx] =Vt[i];
    ut[i+2*nx]=Et[i];
    ut[i+3*nx]=It[i];
    ut[i+4*nx]=Rt[i];
  }
  #
  # Increment calls to flu_1
  ncall <<- ncall+1;
  #
  # Return derivative vector
  return(list(c(ut)));
}

