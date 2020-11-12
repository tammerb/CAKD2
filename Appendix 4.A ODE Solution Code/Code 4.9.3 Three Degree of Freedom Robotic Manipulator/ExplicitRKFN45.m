function [vn,vdn,vddn,Mcond,h,nch]=ExplicitRKFN45(n,tnm,vnm,vdnm,...
    h,nch,par,hvar,nv,Atol,hmax,integ)

%RKFN45 Coefficients
a1=[0,0,0,0,0,0;1/4,0,0,0,0,0;3/32,9/32,0,0,0,0;...
    1932/2197,-7200/2197,7296/2197,0,0,0;439/216,-8,3680/513,-845/4104,0,0;...
    -8/27,2,-3544/2565,1859/4104,-11/40,0];
A=a1*a1;
c=[0;1/4;3/8;12/13;1;1/2];
b1=[25/216,0,1408/2565,2197/4104,-1/5,0];
d1=[16/135,0,6656/12825,28561/56430,-9/50,2/55];
B1=b1*a1;
D1=d1*a1;

%Calculate stage variables

f=ODEfunct(tnm,vnm,vdnm,par,dat,integ);
k1=f;
v2p=vnm+c(2)*h*vdnm;
vd2p=vdnm+h*a1(2,1)*k1;
f=ODEfunct(tnm+c(2)*h,v2p,vd2p,par,dat,integ);
k2=f;
v3p=vnm+c(3)*h*vdnm+(h^2)*A(3,1)*k1;
vd3p=vdnm+h*(a1(3,1)*k1+a1(3,2)*k2);
f=ODEfunct(tnm+c(3)*h,v3p,vd3p,par,dat,integ);
k3=f;
v4p=vnm+c(4)*h*vdnm+(h^2)*(A(4,1)*k1+A(4,2)*k2);
vd4p=vdnm+h*(a1(4,1)*k1+a1(4,2)*k2+a1(4,3)*k3);
f=ODEfunct(tnm+c(4)*h,v4p,vd4p,par,dat,integ);
k4=f;
v5p=vnm+c(5)*h*vdnm+(h^2)*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3);
vd5p=vdnm+h*(a1(5,1)*k1+a1(5,2)*k2+a1(5,3)*k3+a1(5,4)*k4);
f=ODEfunct(tnm+c(5)*h,v5p,vd5p,par,dat,integ);
k5=f;
v6p=vnm+c(6)*h*vdnm+(h^2)*(A(6,1)*k1+A(6,2)*k2+A(6,3)*k3+A(6,4)*k4);
vd6p=vdnm+h*(a1(6,1)*k1+a1(6,2)*k2+a1(6,3)*k3+a1(6,4)*k4+a1(6,5)*k5);
f=ODEfunct(tnm+c(6)*h,v6p,vd6p,par,dat,integ);
k6=f;

k=[k1,k2,k3,k4,k5,k6];

% Evaluate Solution for v, vd, vdd
vn=vnm+h*vdnm+(h^2)*k*B1';
vdn=vdnm+h*k*b1';
f=ODEfunct(tnm+h,vn,vdn,par,dat,integ);
vddn=f;

deriv=1;
[M,gf,M2,gfsv,gfsvd] = AMg(tnm+h,vn,vdn,vddn,par,dat,integ,deriv);
Mcond=cond(M);

if hvar==1  %variable step

%Evaluate imbedded solution for error estimate vnhat and vdnhat
vnhat=vnm+h*vdnm+(h^2)*k*D1';
vdnhat=vdnm+h*k*d1';

%Evaluate  error with sc=Atol, and p=5
diffsqr=(vn-vnhat)'*(vn-vnhat)+(vdn-vdnhat)'*(vdn-vdnhat);
err=sqrt((1/(2*(nv)*(Atol^2)))*diffsqr);

%Change step size

hopt=h*(1/err)^(1/6);

if hopt<h
h=0.5*hopt;
nch=n;
end
if hopt>h
    if n>nch+5
h=2*h;
nch=n;
    end
end

if h>hmax
    h=hmax; 
end
end

end

