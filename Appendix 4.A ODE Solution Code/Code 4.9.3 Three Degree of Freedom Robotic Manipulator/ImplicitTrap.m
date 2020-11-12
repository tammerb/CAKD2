function[vn,vdn,vddn,jiter,R1Norm,JCond]=ImplicitTrap(n,t,...
    V,Vd,Vdd,intol,par,dat,h,integ)

uz=[0;0;1];

% Integration Jacobian Evaluation
vnm=V(:,n-1);
vdnm=Vd(:,n-1);
vddnm=Vdd(:,n-1);
deriv=2;    %Evaluate all derivatives
[M,gf,M2,gfsv,gfsvd] = AMg(t,vnm,vdnm,vddnm,par,dat,integ,deriv);
%M=AM1(vnm,par);
%M2=AM2(vnm,vddnm,par);
%[gsv,gsvd] = Agfsvvd(t,vnm,vdnm,par);
J=M+((h^2)/4)*(M2-gfsv)-(h/2)*gfsvd;
JCond=cond(J);

% Start Integration Step
vdd=vddnm;
if n>2
    vdd=2*vddnm-Vdd(:,n-2);
end

i=1;        %Set iteration counter
err=intol+1;

while err>intol
    
% Residual Calculation
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
vd=vdnm+h*(vddnm+vdd)/2;
deriv = 1;
[M,gf,M2,gfsv,gfsvd] = AMg(t,v,vd,vdd,par,dat,integ,deriv);

R=M*vdd-gf;

if i==1
    R1Norm=norm(R);
end

% Newton Correction
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


