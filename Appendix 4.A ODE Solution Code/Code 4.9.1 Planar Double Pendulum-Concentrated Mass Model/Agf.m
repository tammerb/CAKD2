function g=Agf(t,v,vd,par,dat)

[nv,intol,Atol,hmax,hvar]=BparPart(par);
[m1,m2]=AdatPart(dat);
v1=v(1);
v2=v(2);
vd1=vd(1);
vd2=vd(2);
%Enter right side of ODE, Eq. (4.7.7)

g=[-(m1+m2)*9.8*cos(v1);-m2*9.8*cos(v2)]-...
    [-m2*sin(v2-v1)*(vd2^2);m2*sin(v2-v1)*(vd1^2)];

end

