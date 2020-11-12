function [v,vd,vdd,Lam,ECond,h,nch,jodeiter]=ExplicitRKFN45(n,tnm,...
    vnm,vdnm,vddnm,Lamnm,Uu,U,V,B,q0,h,hmax,par,Atol,nq,nh,nch,hvar)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

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

[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm,vnm,vdnm,vddnm,Lamnm,Uu,...
q0,V,U,B,par);
k1=vdd;
v2p=vnm+c(2)*h*vdnm;
vd2p=vdnm+h*a1(2,1)*k1;
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+c(2)*h,v2p,vd2p,k1,Lam,Uu,...
q0,V,U,B,par);
k2=vdd;
v3p=vnm+c(3)*h*vdnm+(h^2)*A(3,1)*k1;
vd3p=vdnm+h*(a1(3,1)*k1+a1(3,2)*k2);
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+c(3)*h,v3p,vd3p,k2,Lam,Uu,...
q0,V,U,B,par);
k3=vdd;
v4p=vnm+c(4)*h*vdnm+(h^2)*(A(4,1)*k1+A(4,2)*k2);
vd4p=vdnm+h*(a1(4,1)*k1+a1(4,2)*k2+a1(4,3)*k3);
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+c(4)*h,v4p,vd4p,k3,Lam,Uu,...
q0,V,U,B,par);
k4=vdd;
v5p=vnm+c(5)*h*vdnm+(h^2)*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3);
vd5p=vdnm+h*(a1(5,1)*k1+a1(5,2)*k2+a1(5,3)*k3+a1(5,4)*k4);
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+c(5)*h,v5p,vd5p,k4,Lam,Uu,...
q0,V,U,B,par);
k5=vdd;
v6p=vnm+c(6)*h*vdnm+(h^2)*(A(6,1)*k1+A(6,2)*k2+A(6,3)*k3+A(6,4)*k4);
vd6p=vdnm+h*(a1(6,1)*k1+a1(6,2)*k2+a1(6,3)*k3+a1(6,4)*k4+a1(6,5)*k5);
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+c(6)*h,v6p,vd6p,k5,Lam,Uu,...
q0,V,U,B,par);
k6=vdd;

k=[k1,k2,k3,k4,k5,k6];

% Evaluate Solution for v, vd, vdd
v=vnm+h*vdnm+(h^2)*k*B1';
vd=vdnm+h*k*b1';
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+h,v,vd,k6,Lam,Uu,...
q0,V,U,B,par);


if hvar==1  %variable step

%Evaluate imbedded solution for error estimate vnhat and vdnhat
vhat=vnm+h*vdnm+(h^2)*k*D1';
vdhat=vdnm+h*k*d1';

%Evaluate  error with sc=Atol, and p=5
diffsqr=(v-vhat)'*(v-vhat)+(vd-vdhat)'*(vd-vdhat);
err=sqrt((1/(2*(nq-nh)*(Atol^2)))*diffsqr);

%Change step size

hopt=h*(1/err)^(1/6);

if hopt<h
h=hopt;
nch=n;
end
if hopt>h
    if n>nch+5
h=hopt;
nch=n;
    end
end

if h>hmax
    h=hmax; 
end
end

end

