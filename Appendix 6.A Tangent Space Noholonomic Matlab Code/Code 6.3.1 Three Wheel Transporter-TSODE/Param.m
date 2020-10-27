function [q0,qd0,V,W,X,jRepar]=Param(t,q,qd,jRepar,par)
jRepar=jRepar+1;
q0=q;
qd0=qd;

V=eye(3);
C0=Ceval(t,q,par);
W=null(C0);
X=C0';

end

