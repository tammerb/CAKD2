function [qdd,Lam,Qqdd0,LLam0]=Ind3IC(q0,qd0,SMDT,SJDT,STSDAT,par,N)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

%Increment friction coefficients to obtain initial conditions on Lam and
%qdd
S0=SEval(q0,qd0,SMDT,par);
Gam0=GamEval(0,q0,qd0,SJDT,par);
M=MEval(q0,SMDT,par);
Phiq0=PhiqEval(0,q0,SJDT,par);
EE0=[M,Phiq0';Phiq0,zeros(nc)];
EE0Cond=cond(EE0);
Pqdd=[eye(ngc),zeros(ngc,nc)];
PLam=[zeros(nc,ngc),eye(nc)];

%qdd and Lam Estimate with no friction
QA=QAwEval(0,q0,qd0,zeros(nc,1),SMDT,SJDT,STSDAT,par,0,1);
RHS=[QA+S0;-Gam0];
x=EE0\RHS;
qdd=Pqdd*x;
Lam=PLam*x;

w=1;
while w<=N

% Solve for qdd and Lam
i=1;        %Set solution iteration counter
err=1;
while err>intol
%Jacobian Evaluation
QAwsLam=QAwsLamEval(0,q0,qd0,Lam,SJDT,STSDAT,par,w,N);
%[QAsq,QAsqd,QAwsLam]=QAsqqdLam(0,q0,qd0,Lam,SJDT,STSDAT,par);
EE=[M,Phiq0'-QAwsLam;Phiq0,zeros(nc)];
EECond(i)=cond(EE);    
%Residual Calculation
QA=QAwEval(0,q0,qd0,Lam,SMDT,SJDT,STSDAT,par,w,N);
R=[M*qdd+Phiq0'*Lam-QA-S0;Phiq0*qdd+Gam0];
R0norm(i)=norm(R);
% Newton Correction
x=-EE\R;
qdd=qdd+Pqdd*x;
Lam=Lam+PLam*x;

err=norm(R);

i=i+1;

end
Qqdd0(:,w+1)=qdd;
LLam0(:,w+1)=Lam;

w=w+1;
end
end

