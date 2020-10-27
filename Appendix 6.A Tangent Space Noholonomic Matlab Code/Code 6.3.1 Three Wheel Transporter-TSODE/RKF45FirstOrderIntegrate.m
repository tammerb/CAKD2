function [yn,ydn,RM,h,h2,err]=RKF45FirstOrderIntegrate(tnm,ynm,...
    ue,V,W,X,q0,qd0,h,h2,hmax,par)
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);
%RKF45 Coefficients
a=[0,0,0,0,0,0;1/4,0,0,0,0,0;3/32,9/32,0,0,0,0;...
    1932/2197,-7200/2197,7296/2197,0,0,0;439/216,-8,3680/513,-845/4104,0,0;...
    -8/27,2,-3544/2565,1859/4104,-11/40,0];
c=[0;1/4;3/8;12/13;1;1/2];
b1=[25/216,0,1408/2565,2197/4104,-1/5,0];
bh=[16/135,0,6656/12825,28561/56430,-9/50,2/55];

%Calculate stage variables
[f,RM]=ODEfunct(tnm,ynm,ue,V,W,X,q0,qd0,par);
k1=f;
[f,RM]=ODEfunct(tnm+c(2)*h,ynm+h*a(2,1)*k1,...
    ue,V,W,X,q0,qd0,par);
k2=f;
[f,RM]=ODEfunct(tnm+c(3)*h,ynm+h*(a(3,1)*k1+a(3,2)*k2),...
    ue,V,W,X,q0,qd0,par);
k3=f;
[f,RM]=ODEfunct(tnm+c(4)*h,ynm+h*(a(4,1)*k1+a(4,2)*k2+a(4,3)*k3),...
    ue,V,W,X,q0,qd0,par);
k4=f;
[f,RM]=ODEfunct(tnm+c(5)*h,ynm+h*(a(5,1)*k1+a(5,2)*k2+a(5,3)*k3+...
    a(5,4)*k4),ue,V,W,X,q0,qd0,par);
k5=f;
[f,RM]=ODEfunct(tnm+c(6)*h,ynm+h*(a(6,1)*k1+a(6,2)*k2+a(6,3)*k3+...
    a(6,4)*k4+a(6,5)*k5),ue,V,W,X,q0,qd0,par);
k6=f;

%Evaluate yn and ynh and error with Atol=0.001, sc=Atol, and p=5
k=[k1,k2,k3,k4,k5,k6];
yn=ynm+h*k*b1';
ynh=ynm+h*k*bh';


[f,RM]=ODEfunct(tnm+h,yn,ue,V,W,X,q0,qd0,par);
ydn=f;

%Change step size
err=sqrt(((yn-ynh)'*(yn-ynh))/(0.0000000001*(nv+nw)));
hopt=h*(1/err)^(1/6);
h=0.5*hopt;
if h>hmax;
    h=hmax; 
end




end

