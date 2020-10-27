function [v,vd,vdd,q0,U,V,B,jRepar]=Param(n,tnm,Q,Qd,Qdd,SJDT,par,jRepar)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

jRepar=jRepar+1;
q0=Q(:,n-1);
Phiq=PhiqEval(tnm,q0,SJDT,par);
U=Phiq';
B=inv(U'*U);
V=null(Phiq);
v=zeros(nv,1);
qd=Qd(:,n-1);
vd=V'*qd;
qdd=Qdd(:,n-1);
vdd=V'*qdd;

end

