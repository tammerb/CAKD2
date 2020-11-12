function [vd,vdd,q0,U,V,B,jRepar]=Param(n,tnm,Q,Qd,Qdd,PJDT,par,jRepar)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

jRepar=jRepar+1;
q0=Q(:,n-1);
Phiq=PhiqEval(tnm,q0,PJDT,par);
U=Phiq';
B=inv(U'*U);
V=null(U');
qd=Qd(:,n-1);
qdd=Qdd(:,n-1);
vd=V'*qd;
vdd=V'*qdd;

end

