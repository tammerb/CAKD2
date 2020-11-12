function CtLamsq=CtLamsqeval(t,q,Lam,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
Cbarpuz=Cbareval(p,uz);

alp1=Lam(1)*Kbareval(apsa*a,uz)+Lam(2)*Kbareval(bpsa*a,uz)+...
    Lam(3)*eye(4)+Lam(5)*Kbareval(apsa*a,ux)+Lam(6)*Kbareval(apsa*a,uy);
alp2=Lam(1)*Cbarpuz'*apsa+Lam(2)*Cbarpuz'*bpsa+...
    Lam(5)*Cbareval(p,ux)'*apsa+Lam(6)*Cbareval(p,uy)'*apsa;
alp3=Lam(1)*apsa'*Cbarpuz+Lam(2)*bpsa'*Cbarpuz;


CtLamsq=[zeros(3,9);zeros(4,3),alp1,alp2;zeros(2,3),alp3,Lam(4)*eye(2)];


end

