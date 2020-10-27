function R=Resid(v,vd,vdd,Lam,u,par,q0,V,U,B,utol,Btol,h)
%Residual Calculation

[u,Iteru] = usolv(u,v,q0,V,U,B,utol);
q=q0+V*v-U*u;
Phiq = PhiqEval(q);
[B,Biter]=CorrectB(q,B,U,Btol);
D=(eye(11)-U*B*Phiq)*V;
qd=D*vd;
[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
S=SEval(q,qd,par);
QA=QAEval(q,qd,par);
M=MEval(q,par);
R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-S-QA;


end

