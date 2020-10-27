function [v,vd,vdd,Lam,R1n,ImpSoliter,JCond,h,nch]=...
    ImplicitInd0Trap(n,tn,Vv,Vvd,Vvdd,Uu,LLam,q0,V,U,B,...
    h,hmax,nch,PMDT,PTSDAT,PRSDAT,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA,NRSDA]=...
    parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
unm=Uu(:,n-1);
Lamnm=LLam(:,n-1);

Err=2;      %%Criteria for accepting time step
while Err>1

%Jacobian Evaluation
[Rvdd,Rvd,Rv,Phiq,Rhs]=JacobInd0(tnm,vnm,vdnm,vddnm,Lamnm,unm,...
    par,q0,V,U,B,PJDT,PMDT,PTSDAT,PRSDAT);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv,Phiq'];
JCond=cond(J);

%Third derivative calculation, with no explicit time dependence
EE=[Rvdd,Phiq'];
xx=EE\Rhs;
Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];
vdddnm=Pvdd*xx;
Lamdnm=PLam*xx;

vn2=vnm+h*vdnm+((h^2)/2)*vddnm+((h^3)/6)*vdddnm;    %Solution estimates for
vdn2=vdnm+h*vddnm+((h^2)/2)*vdddnm;                 %error control

%Solution Estimates
vdd=vddnm+h*vdddnm;
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
Lam=Lamnm+h*Lamdnm;
u=unm;

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol

%Residual Calculation
R=ResidInd0(tn,v,vd,vdd,Lam,u,B,q0,V,U,PJDT,PMDT,PTSDAT,PRSDAT,par);

if i==1
R1n=norm(R);    
end

x=-J\R;

vdd=vdd+Pvdd*x;
Lam=Lam+PLam*x;

%Evaluate v and vd
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);

err=norm(R);
i=i+1;

end

ImpSoliter=i-1;

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
Er=Er+(vdiff(i)/sciv)^2;        %v-contribution
Er=Er+(vddiff(i)/scivd)^2;      %vd-contribution
i=i+1;
end

Err=sqrt(Er/(2*(ngc-nc)));        %factor of 2 for displ + vel

%Change step size

hopt=h*(1/Err)^(1/3);

if hopt<h
h=0.9*hopt;
nch=n;
if h<10^-5
h=10^-5;
Err=0.5;    
end
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




