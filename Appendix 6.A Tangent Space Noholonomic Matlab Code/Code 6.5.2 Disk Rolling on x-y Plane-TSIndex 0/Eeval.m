function E=Eeval(t,q,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

Bbar=Bbareval(p,apsa*a);

E=[ux',ux'*Bbar,0,0,;uy',uy'*Bbar,0,0];

end

