function S=Seval(q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

p=[q(4);q(5);q(6);q(7)];
pd=[qd(4);qd(5);qd(6);qd(7)];
Gbarpd=Gbareval(pd);

S=[0;0;0;8*Gbarpd'*Jpr*Gbarpd*p;0;0];


end

