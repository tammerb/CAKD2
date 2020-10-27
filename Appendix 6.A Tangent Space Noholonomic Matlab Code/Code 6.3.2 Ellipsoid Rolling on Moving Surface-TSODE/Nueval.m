function Nu=Nueval(t,q,par)
%jEnter expressions for e, Pt, and Nu; example is
%rolling ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
uz=[0;0;1];
gtd=amp*om*sin(om*t);
e=[2*eps*x1;4*eps*y1]*gtd;
Pt=[-gtd*uz;0;0;0;0];
Nu=[-Pt;e];


end

