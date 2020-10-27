function Phi=Phieval(t,q,par,L)
%Evaluate Constraint Function, nhx1; example below is for rolling ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
gt=amp*(1-cos(om*t));
sP=[x1;y1;eps*(x1^2+2*(y1^2))+gt];
w1=[1;0;2*eps*x1];
w2=[0;1;4*eps*y1];
AT=ATran(p);
Phi=[r+AT*apr-sP;w1'*AT*L*apr;w2'*AT*L*apr;...
    0.5*(p'*p-1);0.5*(apr'*L*apr-1)];
end

