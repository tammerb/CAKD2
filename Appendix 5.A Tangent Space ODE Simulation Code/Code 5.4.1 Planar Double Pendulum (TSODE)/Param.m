function [vd,q0,U,V,B,jRepar]=Param(t,n,Q,Qd,jRepar,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

jRepar=jRepar+1;
q0=Q(:,n-1);
Phisq=P1(t,q0,par);
U=Phisq';
B=inv(U'*U);
V=null(U');
qd=Qd(:,n-1);
vd=V'*qd;

end

