function [Rvdd,Rvd,Rv,Phiq,Rhs]=JacobInd0(tn,v,vd,vdd,Lam,u,par,q0,V,U,B,...
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
[QAsq,QAsqd]=QAsqqd(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
P4=P4Eval(tn,q,Lam,PJDT,par);
P21=P2Eval(tn,q,qdd,PJDT,par);
P22=P2Eval(tn,q,qd,PJDT,par);
M2=M2Eval(q,qdd,par);

Rvdd=M*D;
Rvd=-(M*U*B*Gamsqd+QAsqd)*D;
Rv=(M2-M*U*B*P21+P4-M*U*B*Gamsq-QAsq+...
    (M*U*B*Gamsqd+QAsqd)*U*B*(P22+Pstq))*D;

%Right side of third derivative equation 
Rhs=(M*U*B*P21-M2-P4+M*U*B*Gamsq+QAsq)*qd+(M*U*B*Gamsqd+QAsqd)*qdd;

end



