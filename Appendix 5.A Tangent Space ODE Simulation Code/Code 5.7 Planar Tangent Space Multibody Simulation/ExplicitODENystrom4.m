function [v,vd,vdd,ECond]=ExplicitODENystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,npar,PMDT,PTSDAT,PRSDAT,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
u=Uu(:,n-1);

if n-1>npar
u=2*u-Uu(:,n-2);
end

%Nystrom integration step for v and vd plus evaluation of vdd
[vdd,ECond]=ODEfunct(tnm,vnm,vdnm,PMDT,PTSDAT,PRSDAT,PJDT,...
    u,q0,V,U,B,par);
k1=vdd;
[vdd,ECond]=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+(h^2/8)*k1,...
    vdnm+(h/2)*k1,PMDT,PTSDAT,PRSDAT,PJDT,u,q0,V,U,B,par);
k2=vdd;
[vdd,ECond]=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+(h^2/8)*k1,...
    vdnm+(h/2)*k2,PMDT,PTSDAT,PRSDAT,PJDT,u,q0,V,U,B,par);
k3=vdd;
[vdd,ECond]=ODEfunct(tnm+h,vnm+h*vdnm+(h^2/2)*k3,vdnm+h*k3,...
    PMDT,PTSDAT,PRSDAT,PJDT,u,q0,V,U,B,par);
k4=vdd;

% Evaluate Solution for v, vd, vdd
v=vnm+h*vdnm+(h^2/6)*(k1+k2+k3);
vd=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
[vdd,ECond]=ODEfunct(tnm+h,v,vd,PMDT,PTSDAT,PRSDAT,PJDT,...
    u,q0,V,U,B,par);
end

