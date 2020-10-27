function [yn,ydn,RM,h,h2,err]=RKF45FirstOrderIntegrate(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,h,h2,utol,Btol,Htol,par,...
    J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
%RKF45 Coefficients
a=[0,0,0,0,0,0;1/4,0,0,0,0,0;3/32,9/32,0,0,0,0;...
    1932/2197,-7200/2197,7296/2197,0,0,0;439/216,-8,3680/513,-845/4104,0,0;...
    -8/27,2,-3544/2565,1859/4104,-11/40,0];
c=[0;1/4;3/8;12/13;1;1/2];
b1=[25/216,0,1408/2565,2197/4104,-1/5,0];
bh=[16/135,0,6656/12825,28561/56430,-9/50,2/55];

%Calculate stage variables
[f,RM]=ODEfunct(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
k1=f;
[f,RM]=ODEfunct(tnm+c(2)*h,ynm+h*a(2,1)*k1,...
    ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
k2=f;
[f,RM]=ODEfunct(tnm+c(3)*h,ynm+h*(a(3,1)*k1+a(3,2)*k2),...
    ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
k3=f;
[f,RM]=ODEfunct(tnm+c(4)*h,ynm+h*(a(4,1)*k1+a(4,2)*k2+a(4,3)*k3),...
    ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
k4=f;
[f,RM]=ODEfunct(tnm+c(5)*h,ynm+h*(a(5,1)*k1+a(5,2)*k2+a(5,3)*k3+...
    a(5,4)*k4),ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
k5=f;
[f,RM]=ODEfunct(tnm+c(6)*h,ynm+h*(a(6,1)*k1+a(6,2)*k2+a(6,3)*k3+...
    a(6,4)*k4+a(6,5)*k5),ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
k6=f;

%Evaluate yn and ynh and error with Atol=0.001, sc=Atol, and p=5
k=[k1,k2,k3,k4,k5,k6];
yn=ynm+h*k*b1';
ynh=ynm+h*k*bh';


[f,RM]=ODEfunct(tnm+h,yn,ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
ydn=f;

%Change step size
err=sqrt(((yn-ynh)'*(yn-ynh))/(0.0000000001*(nv+nw)));
hopt=h*(1/err)^(1/6);
h=0.5*hopt
if h>hmax
    h=hmax; 
    h2=h^2;
end




end

