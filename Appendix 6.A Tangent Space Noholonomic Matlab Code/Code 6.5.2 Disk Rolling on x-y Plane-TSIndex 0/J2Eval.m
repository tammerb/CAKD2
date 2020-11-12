function J2=J2Eval(t,q,qd,qd0,w,wd,Lam,U,V,W,X,B,H,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

M=Meval(q,par);
[Gamsq,Gamsqd]=Gamsqqdeval(t,q,qd,par);
C=Ceval(t,q,par);
I9=eye(9);
D=(I9-X*H*C)*W;
Gam=Gameval(t,q,qd,par);
x1=-D*wd+X*H*Gam;
C21=C2eval(t,q,x1,par);
Phiq=Phiqeval(t,q,par);
Fbar=(I9-U*B*Phiq)*V;
[Nu,Nusq]=NuNusqeval(t,q,par);
x2=D*w-X*H*C*qd0+qd0+X*H*Nu;
C22=C2eval(t,q,x2,par);
Pbar=-X*H*(C22-Nusq)*Fbar;
CtLamsq=CtLamsqeval(t,q,Lam,par);
[Ssq,Ssqd]=Ssqqdeval(q,qd,par);
[QAsq,QAsqd]=QAsqqdeval(q,qd,par);

R1sv=-V'*Pbar;
R1sw=-V'*D;
R2sv=(M2eval(q,(D*wd-X*H*Gam),par)+M*X*H*C21+CtLamsq-M*X*H*Gamsq-...
    Ssq-QAsq)*Fbar-(M*X*H*Gamsqd+Ssqd+QAsqd)*Pbar;
R2sw=-(M*X*H*Gamsqd+Ssqd+QAsqd)*D;



J2=[R1sv,R1sw,zeros(nv,nh+nd);R2sv,R2sw,zeros(nq,nh+nd)];


end

