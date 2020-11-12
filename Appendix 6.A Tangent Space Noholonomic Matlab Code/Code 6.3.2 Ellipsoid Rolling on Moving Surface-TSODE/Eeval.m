function E=Eeval(t,q,par,L)

[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[r,p,apr,x1,y1]=qPart(q);
w1=[1;0;eps*2*x1];
w2=[0;1;eps*4*y1];
Bbar=BTran(p,apr);
E=[w1',w1'*Bbar,0,0,0,0,0;w2',w2'*Bbar,0,0,0,0,0];

end

