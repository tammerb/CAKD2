function [v,vd,vdd,Lam,ECond,jodeiter1]=ExplicitInd0Nystrom4(n,tn,Vv,Vvd,...
    Vvdd,LLam,Uu,V,U,B,q0,h,npar,SJDT,SMDT,STSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
Lamnm=LLam(:,n-1);
u=Uu(:,n-1);

%Nystrom integration step for v and vd plus evaluation of vdd
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm,vnm,vdnm,vddnm,Lamnm,Uu,...
    q0,V,U,B,SJDT,SMDT,STSDAT,par);
jodeiter1=jodeiter;
k1=vdd;
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm+h/2,vnm+(h/2)*vdnm+(h^2/8)*k1,...
    vdnm+(h/2)*k1,vdd,Lam,Uu,q0,V,U,B,SJDT,SMDT,STSDAT,par);
k2=vdd;
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm+h/2,vnm+(h/2)*vdnm+(h^2/8)*k1,...
    vdnm+(h/2)*k2,vdd,Lam,Uu,q0,V,U,B,SJDT,SMDT,STSDAT,par);
k3=vdd;
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm+h,vnm+h*vdnm+(h^2/2)*k3,...
    vdnm+h*k3,vdd,Lam,Uu,q0,V,U,B,SJDT,SMDT,STSDAT,par);
k4=vdd;

% Evaluate Solution for v, vd, vdd, and Lam
v=vnm+h*vdnm+(h^2/6)*(k1+k2+k3);
vd=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm+h,v,vd,...
    vdd,Lam,Uu,q0,V,U,B,SJDT,SMDT,STSDAT,par);
end






