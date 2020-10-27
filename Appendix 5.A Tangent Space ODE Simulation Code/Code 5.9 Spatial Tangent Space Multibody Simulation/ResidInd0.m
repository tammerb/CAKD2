function R=ResidInd0(tn,v,vd,vdd,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,SJDT,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(tn,q,B,U,SJDT,par);
Phiq=PhiqEval(tn,q,SJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Gam=GamEval(tn,q,qd,SJDT,par);
QA=QAEval(tn,q,qd,SMDT,STSDAT,par);
S=SEval(q,qd,SMDT,par);
M=MEval(q,SMDT,par);

R=M*D*vdd+Phiq'*Lam-(M*U*B*Gam+QA+S);


end

