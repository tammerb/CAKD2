function E3=E3eval(q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[rd,pd,ad]=qdPart(qd);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

Zbarpd=Zbareval(pd,apsa*ad);
Mbarapsa=Mbareval(pd,pd)*apsa;

E3=[0,0,0,ux'*Zbarpd,ux'*Mbarapsa;...
    0,0,0,uy'*Zbarpd,uy'*Mbarapsa];
end

