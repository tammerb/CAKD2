function [q0,qd0,V,U,W,X,B,H,jRepar]=Param(t,q,qd,jRepar,par,L)
jRepar=jRepar+1;
q0=q;
qd0=qd;
Phiq=Phiqeval(t,q,par,L);
V=null(Phiq);
U=Phiq';
B=inv(U'*U);
C0=Ceval(t,q,par,L);
X=C0';
W=null(C0);
H=inv(C0*C0');

end

