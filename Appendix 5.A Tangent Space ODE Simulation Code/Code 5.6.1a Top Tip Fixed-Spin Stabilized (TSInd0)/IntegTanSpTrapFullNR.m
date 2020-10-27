function[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapFullNR(n,npar,Vv,Vvd,Vvdd,...
    LLam,Q,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par)

% Full Jacobian Evaluation
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
Lamnm=LLam(:,n-1);
u=Uu(:,n-1);
qnm=Q(:,n-1);
Phiq = PhiqEval(qnm);

[Rvdd,Rvd,Rv,Phiq]=JacobFull(vnm,vdnm,vddnm,Lamnm,u,par,...
    q0,V,U,B,utol,Btol,h);
J=[Rvdd+(h/2)*Rvd+(h2/4)*Rv,Phiq'];

if n-1==npar
vdd=vddnm;
Lam=Lamnm;
end

if n-1>npar
vdd=2*vddnm-Vvdd(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
u=Uu(:,n-1);

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol;

%ResidualCalculation
R=Resid(v,vd,vdd,Lam,u,par,q0,V,U,B,utol,Btol,h);

if i==1;
    R1n=norm(R);
end

% Full Newton Correction

x=-J\R;

delvdd=[x(1);x(2);x(3)];
delLam=[x(4);x(5);x(6);x(7)];
vdd=vdd+delvdd;
Lam=Lam+delLam;

%Evaluate v and vd
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);

err=norm(R);
i=i+1;

end

jiter=i-1;

end


