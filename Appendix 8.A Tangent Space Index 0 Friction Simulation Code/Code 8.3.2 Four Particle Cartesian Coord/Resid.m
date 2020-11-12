function [R,B,h]=Resid(vdd,Lam,vnm,vdnm,vddnm,u,B,q0,V,U,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);
h=h0;

v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
vd=vdnm+(h/2)*(vddnm+vdd);
[u,uiter]=usolv(u,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par);
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
[Gam,Gamsq,Gamsqd]=GamEval(q,qd,par);

QA=QAEval(q,qd,Lam,par);
M=MEval(q,par);

R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-QA;


end

