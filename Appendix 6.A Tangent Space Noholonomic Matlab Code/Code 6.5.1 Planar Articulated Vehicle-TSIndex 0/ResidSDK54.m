function R=ResidSDK54(t,yk,k,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol)

[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ,h]=Partpar(par);

y=yk+(h/4)*k;
v=[y(1);y(2);y(3);y(4)];
w=y(5);
[u,uiter]=usolv(t,ue,v,q0,V,U,B,utol,par);
q=q0+V*v-U*u;
[B,Biter]=Bcorr(t,q,B,U,Btol,par);
C=Ceval(t,q,par);
[H,Hiter]= Hcorr(H,X,C,Htol,par);
I6=eye(6);
D=(I6-X*H*C)*W;
qd=D*w+(I6-X*H*C)*qd0;
M=Meval(q,par);
Gam=Gameval(t,q,qd,par);
kv=[k(1);k(2);k(3);k(4)];
kw=k(5);

R1=kv-V'*(D*w+(I6-X*H*C)*qd0);
R2=M*D*kw+C'*Lam-M*X*H*Gam;

R=[R1;R2];


end

