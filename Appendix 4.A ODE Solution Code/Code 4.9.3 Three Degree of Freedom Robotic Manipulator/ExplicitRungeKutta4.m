function [vn,vdn,vddn,Mcond]=ExplicitRungeKutta4(tnm,vnm,vdnm,h,par,integ)

%RungeKutta4 integration step for v and vd plus evaluation of vdd
f=ODEfunct(tnm,vnm,vdnm,par,dat,integ);
k1=f;
f=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm,vdnm+(h/2)*k1,par,dat,integ);
k2=f;
f=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/4)*k1,vdnm+(h/2)*k2,par,dat,integ);
k3=f;
f=ODEfunct(tnm+h,vnm+h*vdnm+((h^2)/2)*k2,vdnm+h*k3,par,dat,integ);
k4=f;

% Evaluate Solution for v, vd, vdd
vn=vnm+h*vdnm+((h^2)/6)*(k1+k2+k3);
vdn=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
f=ODEfunct(tnm+h,vn,vdn,par,dat,integ);
vddn=f;
deriv=1;
[M,gf,M2,gfsv,gfsvd] = AMg(tnm+h,vn,vdn,vddn,par,dat,integ,deriv);
Mcond=cond(M);

end

