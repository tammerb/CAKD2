function [R,B,h]=Resid(vdd,Lam,vnm,vdnm,vddnm,u,B,q0,V,U,par,mode,...
    q1b,q2b,q3b)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
vd=vdnm+(h/2)*(vddnm+vdd);
[u,uiter]=usolv(u,v,q0,V,U,B,par,mode,q1b,q2b,q3b);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par,mode);
Phiq=PhiqEval(q,par,mode);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
[Gam,Gamsq,Gamsqd]=GamEval(q,qd,par,mode);

QA=QAEval(q,qd,Lam,par,mode);
M=MEval(q,par);

R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-QA;


end

