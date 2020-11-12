function R=ResidODE(tn,v,vd,vdd,u,B,q0,V,U,PJDT,PMDT,PTSDAT,PRSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Gam=GamEval(tn,q,qd,PJDT,par);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
M=MEval(PMDT,par);

R=D'*M*D*vdd-D'*(M*U*B*Gam+QA);


end

