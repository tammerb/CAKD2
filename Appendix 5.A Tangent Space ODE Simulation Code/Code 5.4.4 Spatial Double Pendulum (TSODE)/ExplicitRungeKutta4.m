function [vn,vdn,vddn,Econd]=ExplicitRungeKutta4(tnm,vnm,vdnm,ue,...
    U,V,B,q0,par,h)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);

%RungeKutta4 integration step for v and vd plus evaluation of vdd
[f,E,B]=ODEfunct(tnm,vnm,vdnm,ue,U,V,B,q0,par);
k1=f;
[f,E,B]=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm,vdnm+(h/2)*k1,...
    ue,U,V,B,q0,par);
k2=f;
[f,E,B]=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/4)*k1,vdnm+(h/2)*k2,...
    ue,U,V,B,q0,par);
k3=f;
[f,E,B]=ODEfunct(tnm+h,vnm+h*vdnm+((h^2)/2)*k2,vdnm+h*k3,...
    ue,U,V,B,q0,par);
k4=f;

% Evaluate Solution for v, vd, vdd
vn=vnm+h*vdnm+((h^2)/6)*(k1+k2+k3);
vdn=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
[f,E,B]=ODEfunct(tnm+h,vn,vdn,ue,U,V,B,q0,par);
vddn=f;
Econd=cond(E);

end

