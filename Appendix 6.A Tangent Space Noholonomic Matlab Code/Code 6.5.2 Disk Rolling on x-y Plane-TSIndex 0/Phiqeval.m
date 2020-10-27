function Phiq=Phiqeval(t,q,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
uz=[0;0;1];
AT=ATran(p);

Phiq=[uz',uz'*Bbareval(p,apsa*a),uz'*AT*apsa;...
    0,0,0,uz'*Bbareval(p,bpsa*a),uz'*AT*bpsa;...
    0,0,0,p',0,0;zeros(1,7),a'];
end
