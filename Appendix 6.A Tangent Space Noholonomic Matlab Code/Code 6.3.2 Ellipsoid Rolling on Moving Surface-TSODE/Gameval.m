function Gam=Gameval(t,q,qd,par,L)

[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
[rd,pd,aprd,x1d,y1d]=qdPart(qd);
gtd=amp*om*sin(om*t);
gtdd=amp*(om^2)*cos(om*t);
uz=[0;0;1];
sPt=gtd*uz;
sPtt=gtdd*uz;
et=eps*[2*x1;4*y1]*gtdd;
eqqd=eps*[2*x1d;4*y1d]*gtd;
Ptt=[-sPtt;0;0;0;0];
Gam=[P2eval(t,q,qd,par,L)+Ptt;E2(q,qd,par)-eqqd-et];



end



