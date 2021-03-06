function [vn,wn,vdn,wdn,Lam,RM]=RungeKutta4FirstOrderIntegrate(tnm,...
    vnm,wnm,ue,V,U,W,X,B,H,q0,qd0,par)
%RungeKutta4 first order integration step for y plus evaluation of yd

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

ynm=[vnm;wnm];

[f,Lam,RM]=ODEfunct(tnm,ynm,ue,V,U,W,X,B,H,q0,qd0,par);
k1=f;
[f,Lam,RM]=ODEfunct(tnm+h/2,ynm+(h/2)*k1,ue,V,U,W,X,B,H,q0,qd0,par);
k2=f;
[f,Lam,RM]=ODEfunct(tnm+h/2,ynm+(h/2)*k2,ue,V,U,W,X,B,H,q0,qd0,par);
k3=f;
[f,Lam,RM]=ODEfunct(tnm+h,ynm+h*k3,ue,V,U,W,X,B,H,q0,qd0,par);
k4=f;

% Evaluate Solution for y

yn=ynm+(h/6)*(k1+2*k2+2*k3+k4);
vn=[yn(1);yn(2);yn(3);yn(4);yn(5)];
wn=[yn(6);yn(7);yn(8)];

[f,Lam,RM]=ODEfunct(tnm+h,yn,ue,V,U,W,X,B,H,q0,qd0,par);
vdn=[f(1);f(2);f(3);f(4);f(5)];
wdn=[f(6);f(7);f(8)];




end

