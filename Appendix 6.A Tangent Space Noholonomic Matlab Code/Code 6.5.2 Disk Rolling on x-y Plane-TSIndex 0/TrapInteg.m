function [vn,wn,vdn,wdn,Lam,RM,jiter,R1n]=TrapInteg(n,t,Vv,Vvd,Ww,Wwd,...
    ue,LLam,V,U,W,X,B,H,q0,qd0,J2,par,nRepar)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

vnm=Vv(:,n-1);
wnm=Ww(:,n-1);
vdnm=Vvd(:,n-1);
wdnm=Wwd(:,n-1);
Lamnm=LLam(:,n-1);

if n-1==nRepar
    vd=vdnm;
    wd=wdnm;
    Lam=Lamnm;    
end

if n-1>nRepar 
    vd=2*vdnm-Vvd(:,n-2);
    wd=2*wdnm-Wwd(:,n-2);
    Lam=2*Lamnm-LLam(:,n-2);
end

vn=vnm+(h/2)*(vdnm+vd);
wn=wnm+(h/2)*(wdnm+wd);

%Jacobian Evaluation

[u,uiter]=usolv(t,ue,vn,q0,V,U,B,par);
q=q0+V*vn-U*u;
M=Meval(q,par);
C=Ceval(t,q,par);
[H,Hiter]= Hcorr(H,X,C,par);
I9=eye(9);
D=(I9-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(t,q,par);
qd=D*wn+(I9-X*H*C)*qd0+X*H*Nu;

J2=J2Eval(t,q,qd,qd0,wn,wd,Lam,U,V,W,X,B,H,par);
Jac=[eye(nv),zeros(nv,nw),zeros(nv,nh+nd);zeros(nq,nv),M*D,C']+(h/2)*J2;
RM=Jac;

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol;

%Residual Calculation
[R,B,H]=ResidTrap(t,vd,wd,Lam,vnm,wnm,vdnm,wdnm,ue,B,H,q0,qd0,...
    V,U,W,X,par);

if i==1;
    R1n=norm(R);
end

% Part Newton Correction

x=-Jac\R;

delvd=[x(1);x(2);x(3);x(4);x(5)];
delwd=[x(6);x(7);x(8)];
delLam=[x(9);x(10);x(11);x(12);x(13);x(14)];
vd=vd+delvd;
wd=wd+delwd;
Lam=Lam+delLam;

%Evaluate w, and v

vn=vnm+(h/2)*(vdnm+vd);
wn=wnm+(h/2)*(wdnm+wd);

err=norm(R);
i=i+1;

end

vdn=vd;
wdn=wd;
jiter=i-1;



end

