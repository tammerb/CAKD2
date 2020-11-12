function g=Agf(t,v,vd,par,dat)

%Enter right side of ODE
[nv,intol,Atol,hmax,hvar]=BparPart(par);
[R,J1,m2,omega0]=AdatPart(dat);

s=sin(v);
c=cos(v);
A=(R*s+(R^2)*s*c/sqrt(4-(R*s)^2));
B=(R*c+(R^2)*(c^2-s^2)/sqrt(4-(R*s)^2)+(R^4)*((s*c)^2)/(4-(R*s)^2)^1.5);

g=-m2*A*B*(vd^2);

end

