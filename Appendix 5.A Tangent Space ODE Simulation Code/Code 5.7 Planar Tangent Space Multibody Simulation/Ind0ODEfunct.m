function [vdd,Lam,ECond]=Ind0ODEfunct(tn,v,vd,PMDT,PTSDAT,PRSDAT,PJDT,...
    u,q0,V,U,B,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA,NRSDA]=...
    parPart(par);

[u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
M=MEval(PMDT,par);
EE=[M*D,Phiq'];
Gam=GamEval(tn,q,qd,PJDT,par);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
RHS=M*U*B*Gam+QA;

x=EE\RHS;

Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];

vdd=Pvdd*x;
Lam=PLam*x;



%ECond=cond(E);
ECond=1;

end





