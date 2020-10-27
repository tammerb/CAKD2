function [v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar)
jRepar=jRepar+1;
q0=Q(:,n-1);
kpr=[0;0;1];
Phiq=PhiqEval(q0);
U=Phiq';
B=inv(U'*U);
V=null(Phiq);
v=[0;0;0];
qd0=Qd(:,n-1);
qdd0=Qdd(:,n-1);
vd=V'*qd0;
vdd=V'*qdd0;

end

