function[v,vd,vdd,Lam,jiter,R1n,h,JCond,Rnorm]=ImplicitTanSpTrap(n,npar,...
    Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);
h=h0;

vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);

if n-1==npar
vdd=vddnm;
Lam=LLam(:,n-1);
end

if n-1>npar
vdd=2*vddnm-Vvdd(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
u=Uu(:,n-1);

%Jacobian Evaluation
u=Uu(:,n-1);
[u,Iteru] = usolv(u,v,q0,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
M=MEval(q,par);
[B,Biter]=CorrectB(q,B,U,par);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par);
[Rvdd,Rvd,Rv,Phiq]=JacobFull(v,vd,vdd,Lam,u,par,q0,V,U,B);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv,Phiq'-QAsLam];
JCond=cond(J);
% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol

    
%Residual Calculation
[R,B,h]=Resid(vdd,Lam,vnm,vdnm,vddnm,u,B,q0,V,U,par);


if i==1
    R1n=norm(R);
end

% Newton Correction

x=-J\R;

delvdd=[x(1);x(2)];
delLam=[x(3);x(4);x(5)];
vdd=vdd+delvdd;
Lam=Lam+delLam;

%Evaluate v and vd
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);

err=norm(R);
i=i+1;


end

jiter=i-1;
Rnorm=norm(R);


end


