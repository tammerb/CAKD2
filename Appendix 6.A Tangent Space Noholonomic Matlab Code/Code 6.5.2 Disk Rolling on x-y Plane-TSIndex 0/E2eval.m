function E2=E2eval(q,x,par)
%Enter terms for (Ex)q
[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[xr,xp,xa]=xPart(x);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

Bbarxp=Bbareval(xp,apsa*a);
Mbarxp=Mbareval(p,xp);

E2=[0,0,0,ux'*Bbarxp,ux'*Mbarxp*apsa;0,0,0,uy'*Bbarxp,uy'*Mbarxp*apsa];

end

