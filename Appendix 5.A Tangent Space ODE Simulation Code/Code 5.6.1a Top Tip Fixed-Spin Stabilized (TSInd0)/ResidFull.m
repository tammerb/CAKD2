function [R,Rvdd,Rvd,Rv,Phiq]=ResidFull(v,vd,vdd,Lam,u,par,...
    q0,V,U,B,utol,Btol,h)
%Residual Calculation

[u,Iteru] = usolv(u,v,q0,V,U,B,utol);
q=q0+V*v-U*u;
Phiq = PhiqEval(q);
[B,Biter]=CorrectB(q,B,U,Btol);
D=(eye(7)-U*B*Phiq)*V;
qd=D*vd;
[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
S=SEval(q,qd,par);
QA=QAEval(q,qd,par);
M=MEval(q,par);
R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-S-QA;


P20=P2Eval(q,qd);
[Sq,Sqd]=SqqdEval(q,qd,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
PqTLsq = PhiqTLamsq(Lam);

Rvdd=M*D;
Rvd=-(M*U*B*Gamsqd+Sqd)*D;
Rv=(M2Eval(q,V*vdd,par)-M*U*B*P2Eval(q,V*vdd)-M2Eval(q,U*B*Gam,par)...
    +M*U*B*P2Eval(q,U*B*Gam)-M*U*B*Gamsq+PqTLsq-Sq...
    +(M*U*B*Gamsqd+Sqd)*U*B*P20-M2Eval(q,U*B*Phiq*V*vdd,par)...
    +M*U*B*P2Eval(q,U*B*Phiq*V*vdd))*D;


end

