function Phiq=Phiqeval(t,q,par,L)
%Enter Jacobian of holonomic constraint functions; example is rolling
%ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
w1=[1;0;2*eps*x1];
w2=[0;1;4*eps*y1];
sPx1=[1;0;2*eps*x1];
sPy1=[0;1;4*eps*y1];
AT=ATran(p);
Bbar=BTran(p,apr);
BbarL=BTran(p,L*apr);
uz=[0;0;1];
Phiq=[eye(3),Bbar,AT,-sPx1,-sPy1;...
    0,0,0,w1'*BbarL,w1'*AT*L,2*eps*apr'*L*AT'*uz,0;...
    0,0,0,w2'*BbarL,w2'*AT*L,0,4*eps*apr'*L*AT'*uz;...
   0,0,0,p',0,0,0,0,0;...
   0,0,0,0,0,0,0,apr'*L,0,0];
end
