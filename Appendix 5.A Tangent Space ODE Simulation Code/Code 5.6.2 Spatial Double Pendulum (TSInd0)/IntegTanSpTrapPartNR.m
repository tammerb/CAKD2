function[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapPartNR(n,npar,...
    Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par,Jvd)

vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);

%Jacobian Evaluation
u=Uu(:,n-1);
[u,Iteru] = usolv(u,vnm,q0,V,U,B,utol);
q=q0+V*vnm-U*u;
M=MEval(q,par);
Phiq = PhiqEval(q);
[B,Biter]=CorrectB(q,B,U,Btol);
D=(eye(11)-U*B*Phiq)*V;

J=[M*D+Jvd,Phiq'];

if n-1==npar
vdd=vddnm;
Lam=LLam(:,n-1);
end

if n-1>npar
vdd=2*vddnm-Vvdd(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+(h2/4)*(vddnm+vdd);
u=Uu(:,n-1);

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol;

%Residual Calculation
R=Resid(v,vd,vdd,Lam,u,par,q0,V,U,B,utol,Btol,h);

if i==1;
    R1n=norm(R);
end

% Part Newton Correction

x=-J\R;

delvdd=[x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8)];
delLam=[x(9);x(10);x(11)];
vdd=vdd+delvdd;
Lam=Lam+delLam;

%Evaluate v and vd
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+(h2/4)*(vddnm+vdd);

err=norm(R);
i=i+1;

end

jiter=i-1;

end


