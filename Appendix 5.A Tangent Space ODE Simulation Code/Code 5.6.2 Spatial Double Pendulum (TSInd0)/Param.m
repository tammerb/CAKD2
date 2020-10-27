function [v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,par,jRepar)
jRepar=jRepar+1;
q0=Q(:,n-1);
kpr=[0;0;1];
Phiq=PhiqEval(q0);
U=Phiq';
B=inv(U'*U);
V=null(Phiq);
v=[0;0;0;0;0;0;0;0];
qd=Qd(:,n-1);
vd=V'*qd;

QA=QAEval(q0,qd,par);
S=SEval(q0,qd,par);
[Gam,Gamsq,Gamsqd] = GamEval(q0,qd);
M=MEval(q0,par);
Coef=[M*V,Phiq'];
RHS=M*U*B*Gam+S+QA;
x=Coef\RHS;
vdd=[x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8)];
Lam=[x(9);x(10);x(11)];

end

