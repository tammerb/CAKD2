function [ErrRvdd,ErrRvd,ErrRv,RvddEst,RvdEst,RvEst] = JacobianCheck(v,...
    vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h,Rvdd,Rvd,Rv)

vdd1=vdd+[h;0;0;0;0;0;0;0];
vdd2=vdd+[0;h;0;0;0;0;0;0];
vdd3=vdd+[0;0;h;0;0;0;0;0];
vdd4=vdd+[0;0;0;h;0;0;0;0];
vdd5=vdd+[0;0;0;0;h;0;0;0];
vdd6=vdd+[0;0;0;0;0;h;0;0];
vdd7=vdd+[0;0;0;0;0;0;h;0];
vdd8=vdd+[0;0;0;0;0;0;0;h];

R=Resid(v,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);

Rvdd1=Resid(v,vd,vdd1,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd2=Resid(v,vd,vdd2,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd3=Resid(v,vd,vdd3,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd4=Resid(v,vd,vdd4,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd5=Resid(v,vd,vdd5,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd6=Resid(v,vd,vdd6,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd7=Resid(v,vd,vdd7,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvdd8=Resid(v,vd,vdd8,Lam,ue,par,q0,V,U,B,utol,Btol,h);

RvddEst=(1/h)*[Rvdd1-R,Rvdd2-R,Rvdd3-R,Rvdd4-R,Rvdd5-R,Rvdd6-R,...
    Rvdd7-R,Rvdd8-R];

vd1=vd+[h;0;0;0;0;0;0;0];
vd2=vd+[0;h;0;0;0;0;0;0];
vd3=vd+[0;0;h;0;0;0;0;0];
vd4=vd+[0;0;0;h;0;0;0;0];
vd5=vd+[0;0;0;0;h;0;0;0];
vd6=vd+[0;0;0;0;0;h;0;0];
vd7=vd+[0;0;0;0;0;0;h;0];
vd8=vd+[0;0;0;0;0;0;0;h];

Rvd1=Resid(v,vd1,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd2=Resid(v,vd2,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd3=Resid(v,vd3,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd4=Resid(v,vd4,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd5=Resid(v,vd5,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd6=Resid(v,vd6,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd7=Resid(v,vd7,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rvd8=Resid(v,vd8,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);

RvdEst=(1/h)*[Rvd1-R,Rvd2-R,Rvd3-R,Rvd4-R,Rvd5-R,Rvd6-R,...
    Rvd7-R,Rvd8-R];

v1=v+[h;0;0;0;0;0;0;0];
v2=v+[0;h;0;0;0;0;0;0];
v3=v+[0;0;h;0;0;0;0;0];
v4=v+[0;0;0;h;0;0;0;0];
v5=v+[0;0;0;0;h;0;0;0];
v6=v+[0;0;0;0;0;h;0;0];
v7=v+[0;0;0;0;0;0;h;0];
v8=v+[0;0;0;0;0;0;0;h];

Rv1=Resid(v1,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv2=Resid(v2,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv3=Resid(v3,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv4=Resid(v4,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv5=Resid(v5,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv6=Resid(v6,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv7=Resid(v7,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);
Rv8=Resid(v8,vd,vdd,Lam,ue,par,q0,V,U,B,utol,Btol,h);

RvEst=(1/h)*[Rv1-R,Rv2-R,Rv3-R,Rv4-R,Rv5-R,Rv6-R,...
    Rv7-R,Rv8-R];


ErrRvdd=Rvdd-RvddEst;
ErrRvd=Rvd-RvdEst;
ErrRv=Rv-RvEst;


end

