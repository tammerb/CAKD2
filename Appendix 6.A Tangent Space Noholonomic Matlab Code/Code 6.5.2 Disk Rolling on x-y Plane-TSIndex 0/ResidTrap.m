function[R,B,H]=ResidTrap(t,vd,wd,Lam,vnm,wnm,vdnm,wdnm,ue,B,H,q0,qd0,...
    V,U,W,X,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

v=vnm+(h/2)*(vdnm+vd);
w=wnm+(h/2)*(wdnm+wd);
[u,uiter]=usolv(t,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Bcorr(t,q,B,U,par);
C=Ceval(t,q,par);
[H,Hiter]= Hcorr(H,X,C,par);
I9=eye(9);
D=(I9-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(t,q,par);
qd=D*w+(I9-X*H*C)*qd0+X*H*Nu;
M=Meval(q,par);
Gam=Gameval(t,q,qd,par);
S=Seval(q,qd,par);
QA=QAEval(t,q,qd,par);

R1=vd-V'*(D*w+(I9-X*H*C)*qd0+X*H*Nu);
R2=M*D*wd+C'*Lam-M*X*H*Gam-S-QA;

R=[R1;R2];


end

