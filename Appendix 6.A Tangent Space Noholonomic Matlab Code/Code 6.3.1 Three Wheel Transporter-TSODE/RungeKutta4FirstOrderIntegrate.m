function [yn,ydn,RM]=RungeKutta4FirstOrderIntegrate(tnm,ynm,...
    ue,V,W,X,q0,qd0,h,h2,par)
%RungeKutta4 first order integration step for y plus evaluation of vdd
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);

[f,RM]=ODEfunct(tnm,ynm,ue,V,W,X,q0,qd0,par);
k1=f;
[f,RM]=ODEfunct(tnm+h/2,ynm+(h/2)*k1,ue,V,W,X,q0,qd0,par);
k2=f;
[f,RM]=ODEfunct(tnm+h/2,ynm+(h/2)*k2,ue,V,W,X,q0,qd0,par);
k3=f;
[f,RM]=ODEfunct(tnm+h,ynm+h*k3,ue,V,W,X,q0,qd0,par);
k4=f;

% Evaluate Solution for y

yn=ynm+(h/6)*(k1+2*k2+2*k3+k4);

[f,RM]=ODEfunct(tnm+h,yn,ue,V,W,X,q0,qd0,par);
ydn=f;



end

