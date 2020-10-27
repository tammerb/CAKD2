function [ErrRvdd,ErrRvd,ErrRv,RvddEst,RvdEst,RvEst] = JacobianCheck(v,...
    vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h,Rvdd,Rvd,Rv)



vdd1=vdd+[h;0;0];
vdd2=vdd+[0;h;0];
vdd3=vdd+[0;0;h];

R=Resid(v,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);

Rvdd1=Resid(v,vd,vdd1,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd2=Resid(v,vd,vdd2,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd3=Resid(v,vd,vdd3,Lam,ue,par,q0,V,U,B,utol,Btol,h);

RvddEst=(1/h)*[Rvdd1-R,Rvdd2-R,Rvdd3-R];

h=h^2;

vd1=vd+[h;0;0];
vd2=vd+[0;h;0];
vd3=vd+[0;0;h];

Rvd1=Resid(v,vd1,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd2=Resid(v,vd2,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd3=Resid(v,vd3,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);

RvdEst=(1/h)*[Rvd1-R,Rvd2-R,Rvd3-R];

h=h^2;

v1=v+[h;0;0];
v2=v+[0;h;0];
v3=v+[0;0;h];

Rv1=Resid(v1,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv2=Resid(v2,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv3=Resid(v3,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);

RvEst=(1/h)*[Rv1-R,Rv2-R,Rv3-R];


ErrRvdd=Rvdd-RvddEst;
ErrRvd=Rvd-RvdEst;
ErrRv=Rv-RvEst;


end

