function [v,vd,vdd,R1n,Jiter,JCond,Jinv,Jinviter,h,nch]=...
    ImplicitODETrap(n,tn,npar,...
    Vv,Vvd,Vvdd,Uu,q0,V,U,B,h,hmax,nch,SMDT,STSDAT,SJDT,par,InvJ,Jinv)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
unm=Uu(:,n-1);

Err=2;      %%Criteria for accepting time step
while Err>1

%Jacobian Evaluation
[Rvdd,Rvd,Rv,Rhs]=Jacob(tnm,vnm,vdnm,vddnm,unm,par,q0,V,U,B,...
    SJDT,SMDT,STSDAT);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv];
JNorm=norm(J);

%Third derivative calculation, with no explicit time dependence
vdddnm=Rvdd\Rhs;

vn2=vnm+h*vdnm+((h^2)/2)*vddnm+((h^3)/6)*vdddnm;    %Solution estimate for 
vdn2=vdnm+h*vddnm+((h^2)/2)*vdddnm;                 %error control                                      


vdd=vddnm+h*vdddnm; %Solution estimate for integration
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
u=unm;

if InvJ==1      %Compute and use Jinv
if n==npar+1
Jinv=inv(J);
Jinviter=0;
JinvNorm=norm(Jinv);
end

if n>npar+1
Jerr=Btol+1;
i=1;
while Jerr>Btol
Jinv=2*Jinv-Jinv*J*Jinv;
Jerr=norm(J*Jinv-eye(nv));
i=i+1;   
end
Jinviter=i-1;
JinvNorm=norm(Jinv);
end
JCond=JinvNorm*JNorm;
end      %End Inveert J

if InvJ==2  %Evaluate Data without Jinv 
JCond=cond(J);
Jinviter=0;
JinvNorm=0;       
end

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol

%Residual Calculation
R=ResidODE(tn,v,vd,vdd,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R1n=norm(R);
end

% Newton Correction

if InvJ==1
delvdd=-Jinv*R;
end

if InvJ==2
delvdd=-J\R;
end

vdd=vdd+delvdd;

%Evaluate v and vd
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);

err=norm(R);
i=i+1;

end
Jiter=i-1;

if hvar==1  %variable step

%Evaluate  error with sciv=Atol*(1+abs(v(i))), scivd=Atol*(1+abs(vd(i)))...
%and p=2
vdiff=v-vn2;
vddiff=vd-vdn2;
Er=0;
i=1;
while i<=ngc-nc
sciv=Atol*(1+abs(v(i)));
scivd=Atol*(1+abs(vd(i)));
Er=Er+(vdiff(i)/sciv)^2+(vddiff(i)/scivd)^2;
i=i+1;
end
Err=sqrt(Er/(2*(ngc-nc)));

%Change step size

hopt=h*(1/Err)^(1/3);

if hopt<h
h=0.9*hopt;
nch=n;
end
if hopt>h
    if n>nch+5
h=min([2*h;0.9*hopt]);
nch=n;
    end
end

if h>hmax
    h=hmax; 
end

end

if hvar==2
Err=0.5;
end

end


end


