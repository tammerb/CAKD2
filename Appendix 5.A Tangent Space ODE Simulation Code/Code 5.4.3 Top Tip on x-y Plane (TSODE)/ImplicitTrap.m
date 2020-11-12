function[vn,vdn,vddn,jiter,R1Norm,JCond,J]=ImplicitTrap(n,tn,...
    Vv,Vvd,Vvdd,Uu,Q,Qd,Qdd,q0,U,V,B,intol,par,h)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

uz=[0;0;1];

% Integration Jacobian Evaluation
q=Q(:,n-1);
qd=Qd(:,n-1);
qdd=Qdd(:,n-1);
[Rv,Rvd,Rvdd]=Jacob(tn,q,qd,qdd,U,V,B,par);
J=Rvdd+h*Rvd/2+(h^2)*Rv/4;
JCond=cond(J);

% Start Integration Step
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);
qddnm=Qdd(:,n-1);
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=V'*qddnm;
vdd=vddnm;

ue=Uu(:,n-1);

if n>2
    vdd=2*vddnm-V'*Qdd(:,n-2);
end

i=1;        %Set iteration counter
err=intol+1;

while err>intol
    
% Residual Calculation
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
vd=vdnm+h*(vddnm+vdd)/2;

R=Resid(tn,ue,v,vd,vdd,q0,V,U,B,par);

if i==1
    R1Norm=norm(R);
end

% Quasi Newton Correction
delvdd=-J\R;
vdd=vdd+delvdd;
err=norm(R);
i=i+1;

end
vddn=vdd;
jiter=i-1;

%Evaluate vn and vdn
vn=vnm+h*vdnm+(h^2)*(vddnm+vdd)/4;
vdn=vdnm+h*(vddnm+vdd)/2;


end


