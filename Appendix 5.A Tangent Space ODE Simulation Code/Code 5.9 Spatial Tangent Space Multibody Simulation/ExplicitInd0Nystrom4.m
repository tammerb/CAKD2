function [v,vd,vdd,Lam,ECond]=ExplicitInd0Nystrom4(n,tn,Vv,Vvd,Uu,...
    V,U,B,q0,h,SMDT,STSDAT,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
u=Uu(:,n-1);

%Nystrom integration step for v and vd plus evaluation of vdd and Lam
[k1,Lam,ECond]=Ind0ODEfunct(tnm,vnm,vdnm,SMDT,STSDAT,SJDT,u,q0,V,U,B,par);
[k2,Lam,ECond]=Ind0ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,...
    vdnm+(h/2)*k1,SMDT,STSDAT,SJDT,u,q0,V,U,B,par);
[k3,Lam,ECond]=Ind0ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,...
    vdnm+(h/2)*k2,SMDT,STSDAT,SJDT,u,q0,V,U,B,par);
[k4,Lam,ECond]=Ind0ODEfunct(tnm+h,vnm+h*vdnm+((h^2)/2)*k3,vdnm+h*k3,...
    SMDT,STSDAT,SJDT,u,q0,V,U,B,par);

% Evaluate Solution for v, vd, vdd, and Lam
v=vnm+h*vdnm+(h^2/6)*(k1+k2+k3);
vd=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
[vdd,Lam,ECond]=Ind0ODEfunct(tnm+h,v,vd,...
    SMDT,STSDAT,SJDT,u,q0,V,U,B,par);
end



