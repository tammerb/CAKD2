function [v,vd,vdd,CondE,jiter]=IntegrateNystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,h2,utol,Btol,par)

jiter=1;

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
u=Uu(:,n-1);

%Nystrom integration step for v and vd plus evaluation of vdd
[f,Lam,B,E]=ODEfunct(tnm,vnm,vdnm,u,V,U,B,q0,utol,Btol,par);
k1=f;
[f,Lam,B,E]=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+(h2/8)*k1,vdnm+(h/2)*k1,...
    u,V,U,B,q0,utol,Btol,par);
k2=f;
[f,Lam,B,E]=ODEfunct(tnm+h/2,vnm+(h/2)*vdnm+(h2/8)*k1,vdnm+(h/2)*k2,...
    u,V,U,B,q0,utol,Btol,par);
k3=f;
[f,Lam,B,E]=ODEfunct(tnm+h,vnm+h*vdnm+(h2/2)*k3,vdnm+h*k3,...
    u,V,U,B,q0,utol,Btol,par);
k4=f;

% Evaluate Solution for v, vd, vdd
v=vnm+h*vdnm+(h2/6)*(k1+k2+k3);
vd=vdnm+(h/6)*(k1+2*k2+2*k3+k4);
[f,Lam,B,E]=ODEfunct(tnm+h,v,vd,u,V,U,B,q0,utol,Btol,par);
vdd=f;
CondE=cond(E);
end

