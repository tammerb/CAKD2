function [vdd,Lam,jodeiter,ECond]=ODEfunct(n,tn,v,vd,vdde,Lame,Uu,...
    q0,V,U,B,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

%Jacobian Evaluation
u=Uu(:,n-1);
[u,Iteru] = usolv(u,v,q0,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
M=MEval(q,par);
[B,Biter]=CorrectB(q,B,U,par);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(tn,q,qd,Lame,par);
E=[M*D,Phiq'-QAsLam];
ECond=cond(E);

% Solve for vdd and Lam

i=1;        %Set solution iteration counter
err=intol+1;
vdd=vdde;
Lam=Lame;
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);

while err>intol
    
%Residual Calculation
QA=QAEval(tn,q,qd,Lam,par);
R=M*D*vdd+Phiq'*Lam-M*U*B*Gam-QA;

% Newton Correction

x=-E\R;

delvdd=[x(1);x(2)];
delLam=x(3);
vdd=vdd+delvdd;
Lam=Lam+delLam;

err=norm(R);
i=i+1;


end

jodeiter=i-1;
ECond=cond(E);


end




