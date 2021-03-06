function R=ResidInd0(tn,v,vd,vdd,Lam,u,B,q0,V,U,...
    PJDT,PMDT,PTSDAT,PRSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Gam=GamEval(tn,q,qd,PJDT,par);
[QA,F]=QAEval(tn,q,qd,Lam,PMDT,PJDT,PTSDAT,PRSDAT,par);
M=MEval(PMDT,par);

R=M*D*vdd+Phiq'*Lam-(M*U*B*Gam+QA);


end


