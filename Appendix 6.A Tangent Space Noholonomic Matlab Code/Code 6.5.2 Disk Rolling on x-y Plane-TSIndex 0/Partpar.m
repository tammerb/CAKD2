function [nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par)

%Add components of par above and below
nq=par(1);
nh=par(2);
nd=par(3);
nv=par(4);
nu=par(5);
nw=par(6);
nx=par(7);
m=par(8);
g=par(9);
Jpr=diag([par(10);par(11);par(12)]);
utol=par(13);
Btol=par(14);
Htol=par(15);
intol=par(16);
h=par(17);

end

