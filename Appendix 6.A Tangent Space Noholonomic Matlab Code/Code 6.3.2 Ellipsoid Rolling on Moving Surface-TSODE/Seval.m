function S=Seval(q,qd,par,Jpr)
%Quadratic velocity term is nqx1; example below is for rolling ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
p=[q(4);q(5);q(6);q(7)];
pd=[qd(4);qd(5);qd(6);qd(7)];
Gd=Geval(pd);
S=[0;0;0;8*Gd'*Jpr*Gd*p;0;0;0;0;0];


end

