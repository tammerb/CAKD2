function [Rvdd,Rvd,Rv,Phiq]=JacobFull(t,v,vd,vdd,Lam,u,par,...
    q0,V,U,B)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

[u,Iteru] = usolv(u,v,q0,V,U,B,par);
q=q0+V*v-U*u;
Phiq = PhiqEval(q,par);
[B,Biter]=CorrectB(q,B,U,par);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
M=MEval(q,par);
P20=P2Eval(q,qd);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);
PqTLsq = PhiqTLamsq(Lam);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(t,q,qd,Lam,par);

Rvdd=M*D;
Rvd=-(M*U*B*Gamsqd+QAsqd)*D;
Rv=(-M*U*B*P2Eval(q,D*vdd-U*B*Gam)-M*U*B*Gamsq+PqTLsq...
    -QAsq+(M*U*B*Gamsqd+QAsqd)*U*B*P20)*D;


end

