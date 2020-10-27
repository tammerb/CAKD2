function J2=J2Eval(t,q,qd,w,wd,Lam,U,V,W,X,B,H,q0,qd0,par)

M=Meval(q,par);
[Gamsq,Gamsqd]=Gamsqqdeval(t,q,qd,par);
C=Ceval(t,q,par);
I6=eye(6);
D=(I6-X*H*C)*W;
Gam=Gameval(t,q,qd,par);
x1=X*H*C*W*wd-W*wd+X*H*Gam;
C21=C2eval(t,q,x1,par);
Phiq=Phiqeval(t,q,par);
Fbar=(I6-U*B*Phiq)*V;
x2=W*w-X*H*C*W*w-X*H*C*qd0+qd0;
C22=C2eval(t,q,x2,par);
Pbar=X*H*C22*Fbar;
CtLamsq=CtLamsqeval(t,q,Lam,par);

R1sw=-V'*D;
R1sv=V'*X*H*C22*Fbar;
R2sw=-M*X*H*Gamsqd*D;
R2sv=(M*X*H*C21+CtLamsq-M*X*H*Gamsq)*Fbar-M*X*H*Gamsqd*Pbar;

J2=[R1sv,R1sw,zeros(4,5);R2sv,R2sw,zeros(6,5)];


end

