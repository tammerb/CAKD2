function [Rvdd,Rvd,Rv,Rhs]=JacobODE(tn,v,vd,vdd,u,par,q0,V,U,B,...
    PJDT,PMDT,PTSDAT,PRSDAT)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(tn,q,PJDT,par);
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Gam=GamEval(tn,q,qd,PJDT,par);
[Gamsq,Gamsqd]=GamsqqdEval(tn,q,qd,PJDT,par);
qdd=D*vdd-U*B*Gam;
M=MEval(PMDT,par);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
P4=P4Eval(tn,q,B'*U'*(M*U*B*Gam+QA),PJDT,par);
P21=P2Eval(tn,q,2*D*vdd-U*B*Gam,PJDT,par);
P22=P2Eval(tn,q,D*vd-U*B*Pst,PJDT,par);
M2=M2Eval(q,D*vdd-U*B*Gam,par);
[QAsq,QAsqd]=QAsqqd(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);

Rvdd=D'*M*D;
Rvd=-D'*(M*U*B*Gamsqd+QAsqd)*D;
Rv=D'*(-M*U*B*P21+P4+M2-M*U*B*Gamsq-QAsq+...
    (M*U*B*Gamsqd+QAsqd)*U*B*(P22+Pstq))*D;

Rhs=D'*(M*U*B*P21-M2-P4+M*U*B*Gamsq+QAsq)*qd+D'*(M*U*B*Gamsqd+QAsqd)*qdd;

end

