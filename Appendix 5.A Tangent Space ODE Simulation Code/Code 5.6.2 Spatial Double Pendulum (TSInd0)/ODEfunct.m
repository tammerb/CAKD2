function [f,Lam,B,E]=ODEfunct(t,v,vd,u,V,U,B,q0,utol,Btol,par)

[u,Iteru] = usolv(u,v,q0,V,U,B,utol);
q=q0+V*v-U*u;
Phiq = PhiqEval(q);
[B,Biter]=CorrectB(q,B,U,Btol);
I=eye(11);
D=(I-U*B*Phiq)*V;
qd=D*vd;
M=MEval(q,par);
QA=QAEval(q,qd,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
S=SEval(q,qd,par);
E=[M*D,Phiq'];
RHS=[M*U*B*Gam+S+QA];
x=E\RHS;
f=[x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8)];
Lam=[x(9);x(10);x(11)];


end

