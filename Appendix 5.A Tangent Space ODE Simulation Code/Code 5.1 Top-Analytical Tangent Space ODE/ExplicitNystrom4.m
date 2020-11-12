function [vn,vdn,vddn,Mcond]=ExplicitNystrom4(tnm,vnm,vdnm,h,par,p0,V,...
    Jpr,m,g,Dcf)

%Nystrom integration step for v and vd, plus evaluation of vdd
f=ODEfunct(tnm,vnm,vdnm,par,p0,V,Jpr,m,g,Dcf);
k1=f;
f=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,vdnm+(h/2)*k1,par,...
    p0,V,Jpr,m,g,Dcf);
k2=f;
f=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,vdnm+(h/2)*k2,par,...
p0,V,Jpr,m,g,Dcf);
k3=f;
f=ODEfunct(tnm+h,vnm+h*vdnm+((h^2)/2)*k3,vdnm+h*k3,par,p0,V,Jpr,m,g,Dcf);
k4=f;

% Evaluate Solution for v, vd, vdd
vn=vnm+h*vdnm+((h^2)/6)*(k1+k2+k3);
vdn=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
f=ODEfunct(tnm+h,vn,vdn,par,p0,V,Jpr,m,g,Dcf);
vddn=f;
Mcond=cond(AM(vn,par,p0,V,Jpr,m,g));

end

