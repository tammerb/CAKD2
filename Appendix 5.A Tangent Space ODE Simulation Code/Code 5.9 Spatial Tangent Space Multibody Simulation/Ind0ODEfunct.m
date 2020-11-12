function [vdd,Lam,ECond]=Ind0ODEfunct(tn,v,vd,SMDT,STSDAT,SJDT,...
    u,q0,V,U,B,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,SJDT,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(tn,q,B,U,SJDT,par);
Phiq=PhiqEval(tn,q,SJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
M=MEval(q,SMDT,par);
E=[M*D,Phiq'];
Gam=GamEval(tn,q,qd,SJDT,par);
QA=QAEval(tn,q,qd,SMDT,STSDAT,par);
S=SEval(q,qd,SMDT,par);
RHS=M*U*B*Gam+QA+S;

x=E\RHS;

Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];
vdd=Pvdd*x;
Lam=PLam*x;

ECond=cond(E);

end





