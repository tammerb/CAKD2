function [q,qd,qdd,ECond,h,nch]=ExplicitRKFN45(n,tn,Q,Qd,h,hmax,par,...
    SMDT,STSDAT,SJDT,nch)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA]=parPart(par);

%RKFN45 Coefficients
a1=[0,0,0,0,0,0;1/4,0,0,0,0,0;3/32,9/32,0,0,0,0;...
    1932/2197,-7200/2197,7296/2197,0,0,0;...
    439/216,-8,3680/513,-845/4104,0,0;...
    -8/27,2,-3544/2565,1859/4104,-11/40,0];
A=a1*a1;
c=[0;1/4;3/8;12/13;1;1/2];
b1=[25/216,0,1408/2565,2197/4104,-1/5,0];
d1=[16/135,0,6656/12825,28561/56430,-9/50,2/55];
B1=b1*a1;
D1=d1*a1;

tnm=tn-h;
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);

Err=2;      %%Criteria for accepting time step
while Err>1

%Calculate stage variables

[qdd,Lam,ECond]=ODEfunct(tnm,qnm,qdnm,SMDT,STSDAT,SJDT,par);
k1=qdd;
q2p=qnm+c(2)*h*qdnm;
qd2p=qdnm+h*a1(2,1)*k1;
[qdd,Lam,ECond]=ODEfunct(tnm+c(2)*h,q2p,qd2p,SMDT,STSDAT,SJDT,par);
k2=qdd;
q3p=qnm+c(3)*h*qdnm+(h^2)*A(3,1)*k1;
qd3p=qdnm+h*(a1(3,1)*k1+a1(3,2)*k2);
[qdd,Lam,ECond]=ODEfunct(tnm+c(3)*h,q3p,qd3p,SMDT,STSDAT,SJDT,par);
k3=qdd;
q4p=qnm+c(4)*h*qdnm+(h^2)*(A(4,1)*k1+A(4,2)*k2);
qd4p=qdnm+h*(a1(4,1)*k1+a1(4,2)*k2+a1(4,3)*k3);
[qdd,Lam,ECond]=ODEfunct(tnm+c(4)*h,q4p,qd4p,SMDT,STSDAT,SJDT,par);
k4=qdd;
q5p=qnm+c(5)*h*qdnm+(h^2)*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3);
qd5p=qdnm+h*(a1(5,1)*k1+a1(5,2)*k2+a1(5,3)*k3+a1(5,4)*k4);
[qdd,Lam,ECond]=ODEfunct(tnm+c(5)*h,q5p,qd5p,SMDT,STSDAT,SJDT,par);
k5=qdd;
q6p=qnm+c(6)*h*qdnm+(h^2)*(A(6,1)*k1+A(6,2)*k2+A(6,3)*k3+A(6,4)*k4);
qd6p=qdnm+h*(a1(6,1)*k1+a1(6,2)*k2+a1(6,3)*k3+a1(6,4)*k4+a1(6,5)*k5);
[qdd,Lam,ECond]=ODEfunct(tnm+c(6)*h,q6p,qd6p,SMDT,STSDAT,SJDT,par);
k6=qdd;

k=[k1,k2,k3,k4,k5,k6];

% Evaluate Solution for v, vd, vdd
q=qnm+h*qdnm+(h^2)*k*B1';
qd=qdnm+h*k*b1';
[qdd,Lam,ECond]=ODEfunct(tnm+h,q,qd,SMDT,STSDAT,SJDT,par);

%Evaluate imbedded solution for error estimate vnhat and vdnhat
qnhat=qnm+h*qdnm+(h^2)*k*D1';
qdnhat=qdnm+h*k*d1';

%Evaluate  error with sciq=Atol*(1+abs(q(i))and sciqd=Atol*(1+abs(qd(i)),
%p=5

qdiff=q-qnhat;
qddiff=qd-qdnhat;

Er=0;
i=1;
while i<=ngc
sciq=Atol*(1+abs(q(i)));
sciqd=Atol*(1+abs(qd(i)));
Er=Er+(qdiff(i)/sciq)^2;    %Position error component
Er=Er+(qddiff(i)/sciqd)^2;  %Velocity error component
i=i+1;
end
Err=sqrt(Er/(2*ngc));    %Factor of 2 removed in pos err cntrl

hopt=h*(1/Err)^(1/6);

%Change step size
if hvar==1  %variable step

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


