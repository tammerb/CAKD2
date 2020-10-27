function [v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);


jRepar=jRepar+1;
q0=Q(:,n-1);
Phiq=PhiqEval(q0,par,mode);
U=Phiq';
B=inv(U'*U);
V=null(Phiq);

if mode==1
v=zeros(nv1,1);
end
if mode>1
v=zeros(nv2,1);
end


qd0=Qd(:,n-1);
qdd0=Qdd(:,n-1);
vd=V'*qd0;
vdd=V'*qdd0;


end

