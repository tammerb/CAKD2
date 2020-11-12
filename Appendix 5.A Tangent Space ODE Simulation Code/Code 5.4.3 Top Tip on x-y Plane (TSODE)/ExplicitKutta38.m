function [vn,vdn,vddn,Econd]=ExplicitKutta38(tnm,vnm,vdnm,ue,...
    U,V,B,q0,par,h)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

%Kutta3/8 integration step for v and vd plus evaluation of vdd
[f,E,B]=ODEfunct(tnm,vnm,vdnm,ue,U,V,B,q0,par);
k1=f;
[f,E,B]=ODEfunct(tnm+h/3,vnm+(h/3)*vdnm,vdnm+(h/3)*k1,...
    ue,U,V,B,q0,par);
k2=f;
[f,E,B]=ODEfunct(tnm+2*h/3,vnm+(2*h/3)*vdnm+(h^2)*k1/3,vdnm+h*(-k1/3+k2),...
    ue,U,V,B,q0,par);
k3=f;
[f,E,B]=ODEfunct(tnm+h,vnm+h*vdnm+(h^2)*(-2*k1/3+k2),vdnm+h*(k1-k2+k3),...
    ue,U,V,B,q0,par);
k4=f;

% Evaluate Solution for v, vd, vdd
vn=vnm+h*vdnm+((h^2)/8)*(k1+2*k2+k3);
vdn=vdnm+(h/8)*(k1+3*k2+3*k3+k4);
[f,E,B]=ODEfunct(tnm+h,vn,vdn,ue,U,V,B,q0,par);
vddn=f;
Econd=cond(E);
end

