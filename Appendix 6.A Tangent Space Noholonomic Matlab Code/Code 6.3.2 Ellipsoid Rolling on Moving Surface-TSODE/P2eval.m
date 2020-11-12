function P2=P2eval(t,q,qd,par,L)
%Evaluate (phiq*x)q*qd
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
[rd,pd,aprd,x1d,y1d]=qdPart(qd);
uz=[0;0;1];
w1=[1;0;2*eps*x1];
w2=[0;1;4*eps*y1];
AT=ATran(p);
Bbard1=BTran(pd,apr);
Bbard2=BTran(pd,aprd);
BbarLd1=BTran(pd,L*apr);
BbarLd2=BTran(pd,L*aprd);
Mbar = MTran(p,pd);
c1=4*eps*x1d*uz'*(BbarLd1*p+AT*L*aprd);
c2=8*eps*y1d*uz'*(BbarLd1*p+AT*L*aprd);
P2=[Bbard1*pd+Bbard2*p+Mbar*aprd-2*eps*(x1d^2)*uz-4*eps*(y1d^2)*uz;...
    w1'*(BbarLd1*pd+BbarLd2*p+Mbar*L*aprd)+c1;...
    w2'*(BbarLd1*pd+BbarLd2*p+Mbar*L*aprd)+c2;pd'*pd;aprd'*L*aprd];
end

