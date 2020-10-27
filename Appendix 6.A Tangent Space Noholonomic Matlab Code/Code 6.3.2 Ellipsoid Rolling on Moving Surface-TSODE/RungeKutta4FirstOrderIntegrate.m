function [yn,ydn,RM]=RungeKutta4FirstOrderIntegrate(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,h,h2,utol,Btol,Htol,par,L,Jpr)
%RungeKutta4 first order integration step for y plus evaluation of vdd
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[f,RM]=ODEfunct(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,utol,Btol,Htol,par,L,Jpr);
k1=f;
[f,RM]=ODEfunct(tnm+h/2,ynm+(h/2)*k1,...
    ue,V,U,W,X,B,H,q0,qd0,utol,Btol,Htol,par,L,Jpr);
k2=f;
[f,RM]=ODEfunct(tnm+h/2,ynm+(h/2)*k2,...
    ue,V,U,W,X,B,H,q0,qd0,utol,Btol,Htol,par,L,Jpr);
k3=f;
[f,RM]=ODEfunct(tnm+h,ynm+h*k3,...
    ue,V,U,W,X,B,H,q0,qd0,utol,Btol,Htol,par,L,Jpr);
k4=f;

% Evaluate Solution for y

yn=ynm+(h/6)*(k1+2*k2+2*k3+k4);

[f,RM]=ODEfunct(tnm+h,yn,...
    ue,V,U,W,X,B,H,q0,qd0,utol,Btol,Htol,par,L,Jpr);
ydn=f;



end

