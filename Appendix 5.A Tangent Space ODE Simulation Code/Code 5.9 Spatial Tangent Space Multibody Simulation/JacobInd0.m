function [Rvdd,Rvd,Rv,Rhs,q,D]=JacobInd0(tn,v,vd,vdd,u,par,q0,V,U,B,...
    SJDT,SMDT,STSDAT)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,SJDT,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(tn,q,SJDT,par);
[B,Biter]=CorrectB(tn,q,B,U,SJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Gam=GamEval(tn,q,qd,SJDT,par);
[Gamsq,Gamsqd]=GamsqqdEval(tn,q,qd,SJDT,par);
qdd=D*vdd-U*B*Gam;
M=MEval(q,SMDT,par);
QA=QAEval(tn,q,qd,SMDT,STSDAT,par);
[QAsq,QAsqd]=QAsqqd(tn,q,qd,SMDT,STSDAT,par);
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
S=SEval(q,qd,SMDT,par);
[Ssq,Ssqd]=SqqdEval(q,qd,SMDT,par);
P4=P4Eval(tn,q,B'*U'*(M*U*B*Gam+QA+S),SJDT,par);
P21=P2Eval(tn,q,D*vdd-U*B*Gam,SJDT,par);
P22=P2Eval(tn,q,D*vd-U*B*Pst,SJDT,par);
M2=M2Eval(q,D*vdd-U*B*Gam,SMDT,par);

Rvdd=M*D;
Rvd=-(M*U*B*Gamsqd+QAsqd)*D;
Rv=(M2-M*U*B*P21+P4-M*U*B*Gamsq-QAsq-Ssq+...
    (M*U*B*Gamsqd+QAsqd+Ssqd)*U*B*(P22+Pstq))*D;

Rhs=(M*U*B*P21-M2-P4+M*U*B*Gamsq+QAsq)*qd+(M*U*B*Gamsqd+QAsqd+Ssqd)*qdd;

end




