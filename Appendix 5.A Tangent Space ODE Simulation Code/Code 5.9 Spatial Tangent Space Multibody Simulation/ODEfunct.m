function [vdd,ECond]=ODEfunct(tn,v,vd,SMDT,STSDAT,SJDT,...
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
E=D'*M*D;
Gam=GamEval(tn,q,qd,SJDT,par);
QA=QAEval(tn,q,qd,SMDT,STSDAT,par);
S=SEval(q,qd,SMDT,par);
RHS=D'*(M*U*B*Gam+QA+S);

vdd=E\RHS;

ECond=cond(E);

end




