function P2=P2eval(t,q,x,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[xr,xp,xa]=xPart(x);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

Mbarxp=Mbareval(p,xp);

P2=[0,0,0,uz'*(Bbareval(xp,apsa*a)+Bbareval(p,apsa*xa)),uz'*Mbarxp*apsa;...
    0,0,0,uz'*(Bbareval(xp,bpsa*a)+Bbareval(p,bpsa*xa)),uz'*Mbarxp*bpsa;
    0,0,0,xp',0,0;zeros(1,7),xa'];
end

