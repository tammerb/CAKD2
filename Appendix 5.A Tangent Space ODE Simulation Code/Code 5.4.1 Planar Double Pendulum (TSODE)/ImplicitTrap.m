function[vn,vdn,vddn,jiter,R1Norm,JCond]=ImplicitTrap(n,tn,...
    Vv,Vvd,Vvdd,Uu,Q,Qd,Qdd,q0,U,V,B,intol,par,h)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

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

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);
R=D'*M*D*vdd-D'*(M*U*B*Gamma+S+QA);

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


