function [q,qd,qdd,qdddnm,Lam,Lamdnm,R1n,Jiter,JCond,h,nch,Err]=...
    ImplicitIndex2(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
    PMDT,PTSDAT,PRSDAT,PJDT,par,alpha,nch)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA,NRSDA]=parPart(par);

beta=((1-alpha)^2)/4;
delta=(1-2*alpha)/2;

tnm=tn-h;
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);
qddnm=Qdd(:,n-1);
Lamnm=LLam(:,n-1);

Err=2;      %%Criteria for accepting time step
while Err>1
    
Phiq=PhiqEval(tnm,qnm,PJDT,par);
M=MEval(PMDT,par);
M2=M2Eval(qnm,qddnm,par);
QA=QAEval(tnm,qnm,qdnm,PMDT,PTSDAT,PRSDAT,par);
[QAsq,QAsqd]=QAsqqd(tnm,qnm,qdnm,PMDT,PTSDAT,PRSDAT,par);
P4=P4Eval(tnm,qnm,Lamnm,PJDT,par);
P2qdd=P2Eval(tnm,qnm,qddnm,PJDT,par);
P2qd=P2Eval(tnm,qnm,qdnm,PJDT,par);
[Gamsq,Gamsqd]=GamsqqdEval(tnm,qnm,qdnm,PJDT,par);

%Third Derivative calculation
E=[M,Phiq';Phiq,zeros(nc,nc)];
Rhs=[(-M2-P4+QAsq)*qdnm+QAsqd*qddnm;...
    (-P2qdd-Gamsq)*qdnm-Gamsqd*qddnm];
x=E\Rhs;
Pqdd=[eye(ngc),zeros(ngc,nc)];
PLam=[zeros(nc,ngc),eye(nc)];
qdddnm=Pqdd*x;
Lamdnm=PLam*x;

%Estimate for qdd and Lam
qdd=qddnm+h*qdddnm;
Lam=Lamnm+h*Lamdnm;

qd=qdnm+h*((1-delta)*qddnm+delta*qdd);
q=qnm+h*qdnm+((h^2)/2)*((1-2*beta)*qddnm+2*beta*qdd);

Phiq=PhiqEval(tn,q,PJDT,par);
M=MEval(PMDT,par);
M2=M2Eval(q,qdd,par);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
[QAsq,QAsqd]=QAsqqd(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
P4=P4Eval(tn,q,Lam,PJDT,par);
P2=P2Eval(tn,q,qdd,PJDT,par);
[Gamsq,Gamsqd]=GamsqqdEval(tn,q,qd,PJDT,par);

%Jacobian Evaluation
R1sqdd=(1/(1+alpha))*M-h*delta*QAsqd+...
    ((h^2)*beta)*((1/(1+alpha))*M2+P4-QAsq);
R2sqdd=Phiq+(beta*h/delta)*P2qd;
J=[R1sqdd,Phiq';R2sqdd,zeros(nc,nc)];
JCond=cond(J);

% Solve Discretized Equations
i=1;        %Set solution iteration counter
err=intol+1;
R1pnm=-(alpha/(1+alpha))*(PhiqEval(tnm,qnm,PJDT,par)'*Lamnm-QA);
while err>intol
    
%Residual Calculation
Phiq=PhiqEval(tn,q,PJDT,par);
M=MEval(PMDT,par);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par); 
R=[(1/(1+alpha))*M*qdd+Phiq'*Lam-QA+R1pnm;(1/(delta*h))*Phiq*qd];
if i==1
R1n=norm(R);
end

% Newton Correction
x=-J\R;
qdd=qdd+Pqdd*x;
Lam=Lam+PLam*x;

%Evaluate q and qd
qd=qdnm+h*((1-delta)*qddnm+delta*qdd);
q=qnm+h*qdnm+((h^2)/2)*((1-2*beta)*qddnm+2*beta*qdd);

err=norm(R);
i=i+1;
end
Jiter=i-1;

%Step size control
if hvar==1
%Higher order estimate of solution
qe=qnm+h*qdnm+0.5*(h^2)*qddnm+((h^3)/6)*qdddnm;
qde=qdnm+h*qddnm+0.5*(h^2)*qdddnm;

%Evaluate  error with sci=Atol(1+abs(q(i))+abs(qd(i))), and p=2
qdiff=q-qe;
qddiff=qd-qde;
Er=0;
i=1;
while i<=ngc
sciq=Atol*(1+abs(q(i)));
sciqd=Atol*(1+abs(qd(i)));
Er=Er+(qdiff(i)/sciq)^2;    %Position error component
Er=Er+(qddiff(i)/sciqd)^2;  %Velocity error component
i=i+1;
end
Err=sqrt(Er/(2*ngc));

%Change step size

hopt=h*(1/Err)^(1/3);

if hopt<h
h=0.9*hopt;
nch=n;
end
if hopt>h
    if n>nch+5
h=min([2*h;0.9*hopt]);
nch=n;
    end
end

if h>hmax
    h=hmax; 
end

if h<10^-5
h=10^-5;
Err=0.5;    
end

end

if hvar==2
Err=0.5;
end

end

end





