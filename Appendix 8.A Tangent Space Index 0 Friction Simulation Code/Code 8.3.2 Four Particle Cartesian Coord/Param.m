function [v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

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
vdd=V'*qdd0;

end

