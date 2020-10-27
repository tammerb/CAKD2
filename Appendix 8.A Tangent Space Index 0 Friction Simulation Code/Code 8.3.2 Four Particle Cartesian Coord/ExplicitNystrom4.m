function [v,vd,vdd,Lam,ECond,jodeiter]=ExplicitNystrom4(n,tnm,...
    Vv,Vvd,Vvdd,LLam,Uu,U,V,B,q0,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);
h=h0;

vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vdde=Vvdd(:,n-1);
Lame=LLam(:,n-1);

%Nystrom integration step for v and vd, plus evaluation of vdd
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm,vnm,vdnm,vdde,Lame,Uu,...
    q0,V,U,B,par);
k1=vdd;

[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,...
    vdnm+(h/2)*k1,vdd,Lam,Uu,...
    q0,V,U,B,par);
k2=vdd;

[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+h/2,vnm+(h/2)*vdnm+((h^2)/8)*k1,...
    vdnm+(h/2)*k2,vdd,Lam,Uu,q0,V,U,B,par);
k3=vdd;

[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+h,vnm+h*vdnm+((h^2)/2)*k3,...
    vdnm+h*k3,vdd,Lam,Uu,...
    q0,V,U,B,par);
k4=vdd;

% Evaluate Solution for v, vd, vdd
v=vnm+h*vdnm+((h^2)/6)*(k1+k2+k3);
vd=vdnm+(h/6)*(k1+2*k2+2*k3+k4);

[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tnm+h,v,vd,vdd,Lam,Uu,...
    q0,V,U,B,par);


end

