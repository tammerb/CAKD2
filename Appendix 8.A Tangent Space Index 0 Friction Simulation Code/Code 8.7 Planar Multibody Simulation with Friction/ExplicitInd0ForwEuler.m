function [v,vd,vdd,Lam,ECond,jodeiter]=ExplicitInd0ForwEuler(n,tn,Vv,...
    Vvd,Vvdd,LLam,Uu,V,U,B,q0,h,npar,PJDT,PMDT,PTSDAT,PRSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
Lamnm=LLam(:,n-1);
u=Uu(:,n-1);

%Forward Euler integration step for v and vd plus evaluation of vdd & Lam
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm,vnm,vdnm,vddnm,Lamnm,Uu,...
    q0,V,U,B,PJDT,PMDT,PTSDAT,PRSDAT,par);

v=vnm+h*vdnm+(h^2)*vdd;
vd=vdnm+h*vdd;
[vdd,Lam,jodeiter,ECond]=FrODEfunct(n,tnm+h,v,vd,...
    vdd,Lam,Uu,q0,V,U,B,PJDT,PMDT,PTSDAT,PRSDAT,par);
end



