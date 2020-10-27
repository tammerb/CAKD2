function[R,B,H]=Resid(t,vd,wd,Lam,vnm,wnm,vdnm,wdnm,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol)

[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ,h]=Partpar(par);
v=vnm+(h/2)*(vdnm+vd);
w=wnm+(h/2)*(wdnm+wd);
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

R1=vd-V'*(D*w+(I6-X*H*C)*qd0);
R2=M*D*wd+C'*Lam-M*X*H*Gam;

R=[R1;R2];


end

