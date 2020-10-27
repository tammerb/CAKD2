function [Rvdd,Rvd,Rv,Phiq]=JacobFull(v,vd,vdd,Lam,u,par,...
    q0,V,U,B,mode,q1b,q2b,q3b,Sw21)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

[u,Iteru] = usolv(u,v,q0,V,U,B,par,mode,q1b,q2b,q3b);
q=q0+V*v-U*u;
Phiq = PhiqEval(q,par,mode);
[B,Biter]=CorrectB(q,B,U,par,mode);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
M=MEval(q,par);
P20=P2Eval(q,qd,par,mode);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par,mode);
PqTLsq = PhiqTLamsq(Lam,par,mode);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par,mode);

Rvdd=M*D;
Rvd=-(M*U*B*Gamsqd+QAsqd)*D;
Rv=(-M*U*B*P2Eval(q,D*vdd-U*B*Gam,par,mode)-M*U*B*Gamsq+PqTLsq...
    -QAsq+(M*U*B*Gamsqd+QAsqd)*U*B*P20)*D;


end

