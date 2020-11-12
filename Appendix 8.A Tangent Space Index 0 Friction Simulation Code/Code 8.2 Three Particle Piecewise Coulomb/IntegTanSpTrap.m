function[v1,v11d,v1dd,Lam1,jiter,R1n,h,JCond,Rnorm]=IntegTanSpTrap1(n,npar,...
    Vv1,Vv1d,Vv1dd,LLam1,Uu1,q0,V1,U1,B1,par)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h,mode]=Partpar(par);

v1nm=Vv1(:,n-1);
v1dnm=Vv1d(:,n-1);
v1ddnm=Vv1dd(:,n-1);

if n-1==npar
v1dd=v1ddnm;
Lam1=LLam1(:,n-1);
end

if n-1>npar
v1dd=2*v1ddnm-Vv1dd(:,n-2);
Lam1=2*LLam1(:,n-1)-LLam1(:,n-2);
end

v1d=v1dnm+(h/2)*(v1ddnm+v1dd);
v1=v1nm+h*v1dnm+((h^2)/4)*(v1ddnm+v1dd);
u1=Uu1(:,n-1);

%Jacobian Evaluation
[u1,Iteru] = usolv(u1,v1,q0,V1,U1,B1,par);
q=q0+V1*v1-U1*u1;
Phiq=PhiqEval(q,par);
D1=(eye(nq)-U1*B1*Phiq)*V1;
qd=D1*vd;
M=MEval(q,par);
[B1,Biter]=CorrectB(q,B1,U1,par);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam1,par);
[Rvdd,Rvd,Rv,Phiq]=JacobFull(v1,v1d,v1dd,Lam1,u1,par,q0,V1,U1,B1);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv,Phiq'-QAsLam];
JCond=cond(J);
% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol;

    
%Residual Calculation
[R,B,h]=Resid(v1dd,Lam1,v1nm,v1dnm,v1ddnm,u1,B1,q0,V1,U1,par);


if i==1;
    R1n=norm(R);
end

% Newton Correction

x=-J\R;

delv1dd=[x(1);x(2)];
delLam1=x(3);
v1dd=v1dd+delv1dd;
Lam1=Lam1+delLam1;

%Evaluate v and vd
v1d=v1dnm+(h/2)*(v1ddnm+v1dd);
v1=v1nm+h*v1dnm+((h^2)/4)*(v1ddnm+v1dd);

err=norm(R);
i=i+1;

    if i>6
err=intol/2;
    end

end

jiter=i-1;
Rnorm=norm(R);


end


