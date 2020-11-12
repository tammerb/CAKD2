function [qdd,Lam,Qdd0,LLam0,w]=Ind0IC(q0,qd0,SMDT,SJDT,STSDAT,par,N)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

%Evaluate Initial Acceleration and Lagrange Multipliers;

%Initial Parameterization
Phiq=PhiqEval(0,q0,SJDT,par);
U=Phiq';
B=inv(U'*U);
V=null(U');
v=zeros(nv,1);
vd=V'*qd0;
u=zeros(nu,1);
Vv(:,1)=v;
Vvd(:,1)=vd;
Uu(:,1)=u;

%Increment friction coefficients to obtain initial conditions on Lam and
%vdd
w=0;
Lam=zeros(nc,1);
vdd=zeros(nv,1);
while w<=N
[vdd,Lam,jodeiter,ECond]=FrODEfunct0w(q0,qd0,vdd,Lam,...
   V,U,B,SJDT,SMDT,STSDAT,par,w,N);

w=w+1;

Gam=GamEval(0,q0,qd0,SJDT,par);
qdd=V*vdd-U*B*Gam;
Econd0(w)=ECond;
jodeiter0(w)=jodeiter;
Qdd0(:,w)=qdd;
LLam0(:,w)=Lam;

end

end

