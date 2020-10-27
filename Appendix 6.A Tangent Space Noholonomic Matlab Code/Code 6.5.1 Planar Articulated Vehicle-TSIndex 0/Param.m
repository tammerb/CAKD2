function [q0,qd0,V,U,W,X,B,H,jRepar]=Param(t,q,qd,jRepar,par)
jRepar=jRepar+1;
q0=q;
qd0=qd;
Phiq=Phiqeval(t,q,par);
V=null(Phiq);
U=Phiq';
B=inv(U'*U);
C0=Ceval(t,q,par);
W=null(C0);
X=C0';
H=inv(C0*C0');

end

