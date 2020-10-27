function Phi=Phieval(t,q,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
uz=[0;0;1];

AT=ATran(p);
Phi=[uz'*(r+AT*apsa*a);uz'*AT*bpsa*a;(p'*p-1)/2;(a'*a-1)/2];



end