function [vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tn,v,vd,vdd,Lam,u,...
    q0,V,U,B,SJDT,SMDT,STSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

%Evaluate terms that do not depend on vdd and Lam
[u,Iteru]=usolv(tn,u,v,q0,SJDT,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(tn,q,SJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
qd=D*vd;
M=MEval(q,SMDT,par);
Gam=GamEval(tn,q,qd,SJDT,par);
[B,Biter]=CorrectB(tn,q,B,U,SJDT,par);
Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];

% Solve for vdd and Lam

i=1;        %Set solution iteration counter
err=intol+1;


while err>intol
%Jacobian Evaluation
QAsLam=QAsLamEval(tn,q,qd,Lam,SJDT,STSDAT,par);
E=[M*D,Phiq'-QAsLam];
    
%Residual Calculation
QA=QAEval(tn,q,qd,Lam,SMDT,SJDT,STSDAT,par);
S=SEval(q,qd,SMDT,par);
R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-QA-S;

% Newton Correction
x=-E\R;
vdd=vdd+Pvdd*x;
Lam=Lam+PLam*x;

err=norm(R)+norm(x);
i=i+1;

end

jodeiter=i-1;
ECond=cond(E);

end


