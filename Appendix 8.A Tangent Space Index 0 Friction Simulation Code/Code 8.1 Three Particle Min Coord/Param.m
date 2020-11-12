function [v,vd,vdd,Lam,q0,U,V,B,jRepar]=Param(n,tn,Q,Qd,Qdd,LLam,Uu,...
    par,jRepar)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);
jRepar=jRepar+1;
q0=Q(:,n-1);
Phiq=PhiqEval(q0,par);
U=Phiq';
B=inv(U'*U);
V=null(Phiq);
v=zeros(nv,1);
qd0=Qd(:,n-1);
qdd0=Qdd(:,n-1);
vd=V'*qd0;
vdde=V'*qdd0;
Lame=LLam(:,n-1);
[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tn,v,vd,vdde,Lame,Uu,q0,V,U,B,par);

end

