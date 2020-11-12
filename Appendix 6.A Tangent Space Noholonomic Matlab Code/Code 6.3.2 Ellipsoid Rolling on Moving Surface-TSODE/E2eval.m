function E2=E2eval(q,qd,par)
%Enter terms for (Eqd)q*qd
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
[rd,pd,aprd,x1d,y1d]=qdPart(qd);
w1=[1;0;eps*2*x1];
w2=[0;1;eps*4*y1];
uz=[0;0;1];
Bbar=BTran(p,apr);
Bbard=BTran(pd,apr);
Mbar=MTran(p,pd);
E2=[2*eps*x1d*uz'*(rd+Bbar*pd)+w1'*(Bbard*pd+Mbar*aprd);...
    4*eps*y1d*uz'*(rd+Bbar*pd)+w2'*(Bbard*pd+Mbar*aprd)];

end

