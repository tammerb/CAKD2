function [Rv,Rvd,Rvdd]=JacobEval(n,q,qd,Qd,vdd,Lam,U,V,B,par,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

Phiq=PhiqEval(q,par,mode);
P20=P2Eval(q,qd,par,mode);
QA=QAEval(q,qd,Lam,par,mode);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd]=GamEval(q,qd,par,mode);
PqTLsq = PhiqTLamsq(Lam,par,mode);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par,mode);

Rvdd=M*V;
Rvd=-(M*U*B*Gamsqd+QAsqd)*V;
Rv=(-M*U*B*P2Eval(q,V*vdd-U*B*Gam,par,mode)-M*U*B*Gamsq+PqTLsq-QAsq...
    +(M*U*B*Gamsqd+QAsqd)*U*B*P20)*V;


end



