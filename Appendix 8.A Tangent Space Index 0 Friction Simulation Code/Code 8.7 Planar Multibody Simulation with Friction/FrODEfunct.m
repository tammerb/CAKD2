function [vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tn,v,vd,vdd,Lam,u,...
    q0,V,U,B,PJDT,PMDT,PTSDAT,PRSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

%Evaluate terms that do not depend on vdd and Lam
[u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(tn,q,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
qd=D*vd;
M=MEval(PMDT,par);
Gam=GamEval(tn,q,qd,PJDT,par);
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];

% Solve for vdd and Lam

i=1;        %Set solution iteration counter
err=intol+1;


while err>intol
%Jacobian Evaluation
QAsLam=QAsLamEval(tn,q,qd,Lam,PJDT,par);
E=[M*D,Phiq'-QAsLam];
    
%Residual Calculation
[QA,F]=QAEval(tn,q,qd,Lam,PMDT,PJDT,PTSDAT,PRSDAT,par);
R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-QA;

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

