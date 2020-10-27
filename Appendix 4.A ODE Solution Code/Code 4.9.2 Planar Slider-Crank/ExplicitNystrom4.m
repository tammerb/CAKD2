function [vn,vdn,vddn,Mcond]=ExplicitNystrom4(tnm,vnm,vdnm,h,par,dat)

%Nystrom integration step for v and vd, plus evaluation of vdd
f=ODEfunct(tnm,vnm,vdnm,par,dat);
k1=f;
f=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,vdnm+(h/2)*k1,par,dat);
k2=f;
f=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,vdnm+(h/2)*k2,par,dat);
k3=f;
f=ODEfunct(tnm+h,vnm+h*vdnm+((h^2)/2)*k3,vdnm+h*k3,par,dat);
k4=f;

% Evaluate Solution for v, vd, vdd
vn=vnm+h*vdnm+((h^2)/6)*(k1+k2+k3);
vdn=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
f=ODEfunct(tnm+h,vn,vdn,par,dat);
vddn=f;
Mcond=cond(AM(vn,par,dat));

end
