function P3=P3eval(t,q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[rd,pd,ad]=qdPart(qd);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
uz=[0;0;1];
Mbarpdpd=Mbareval(pd,pd);

P3=[0,0,0,uz'*(Bbareval(pd,apsa*ad)+Zbareval(pd,apsa*ad)),...
    uz'*Mbarpdpd*apsa;
    0,0,0,uz'*(Bbareval(pd,bpsa*ad)+Zbareval(pd,bpsa*ad)),...
    uz'*Mbarpdpd*bpsa;zeros(2,9)];

end

