function [gsv,gsvd] = Agfsvvd(t,v,vd,par,dat)

[nv,intol,Atol,hmax,hvar]=BparPart(par);
[m1,m2]=AdatPart(dat);
v1=v(1);
v2=v(2);
v1d=vd(1);
v2d=vd(2);

%Enter derivatives of g with respect to v and vd
s1=sin(v1);
s2=sin(v2);
c2m1=cos(v2-v1);
s2m1=sin(v2-v1);

gsv=[(m1+m2)*9.8*s1-m2*c2m1*v2d^2,m2*c2m1*v2d^2;...
    m2*c2m1*v1d^2,m2*9.8*s2-m2*c2m1*v1d^2];
gsvd=[0,2*m2*s2m1*v2d;-2*m2*s2m1*v1d,0];

end

