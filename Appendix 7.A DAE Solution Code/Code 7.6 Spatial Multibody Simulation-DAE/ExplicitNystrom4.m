function [q,qd,qdd,Lam,ECond]=ExplicitNystrom4(n,tn,Q,Qd,h,...
    SMDT,STSDAT,SJDT,par)

tnm=tn-h;
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);

%Nystrom integration step for q and qd plus evaluation of qdd
[qdd,Lam,ECond]=ODEfunct(tnm,qnm,qdnm,SMDT,STSDAT,SJDT,par);
k1=qdd;
[qdd,Lam,ECond]=ODEfunct(tnm+h/2,qnm+(h/2)*qdnm+(h^2/8)*k1,...
    qdnm+(h/2)*k1,SMDT,STSDAT,SJDT,par);
k2=qdd;
[qdd,Lam,ECond]=ODEfunct(tnm+h/2,qnm+(h/2)*qdnm+(h^2/8)*k1,...
    qdnm+(h/2)*k2,SMDT,STSDAT,SJDT,par);
k3=qdd;
[qdd,Lam,ECond]=ODEfunct(tnm+h,qnm+h*qdnm+(h^2/2)*k3,...
    qdnm+h*k3,SMDT,STSDAT,SJDT,par);
k4=qdd;

% Evaluate Solution for q, qd, qdd, and Lam
q=qnm+h*qdnm+(h^2/6)*(k1+k2+k3);
qd=qdnm+(h/6)*(k1+2*k2+2*k3+k4);

[qdd,Lam,ECond]=ODEfunct(tnm+h,q,qd,SMDT,STSDAT,SJDT,par);
end

