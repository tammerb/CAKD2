function [vdd,Lam,jodeiter,ECond]=FrODEfunct0w(q0,qd0,vdd,Lam,...
    V,U,B,PJDT,PMDT,PTSDAT,PRSDAT,par,w,N)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

%Evaluate terms that do not depend on vdd and Lam
M=MEval(PMDT,par);
Gam=GamEval(0,q0,qd0,PJDT,par);
Phiq=PhiqEval(0,q0,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];

% Solve for vdd and Lam

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol
%Jacobian Evaluation
QAwsLam=QAwsLamEval(q0,qd0,Lam,PJDT,par,w,N);
E=[M*D,Phiq'-QAwsLam];
    
%Residual Calculation
QAw=QAwEval(q0,qd0,Lam,PMDT,PJDT,PTSDAT,PRSDAT,par,w,N);
R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-QAw;

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


