function [gsv,gsvd] = Agfsvvd(t,v,vd,par)

%Enter derivatives of g with respect to v and vd

[nv,intol,Atol,hmax,hvar]=BparPart(par);
[R,J1,m2,omega0]=AdatPart(dat);

s=sin(v);
c=cos(v);
A=(R*s+(R^2)*s*c/sqrt(4-(R*s)^2));
B=(R*c+(R^2)*(c^2-s^2)/sqrt(4-(R*s)^2)+(R^4)*((s*c)^2)/(4-(R*s)^2)^1.5);
C=(-R*s-2*(R^2)*s*c/sqrt(4-(R*s)^2)+...
    3*(R^4)*(s*c^3-c*s^3)/((4-(R*s)^2))^1.5+...
    3*(R^6)*((s*c)^3)/(2*((4-(R*s)^2))^2.5));

gsv=m2*(B^2+A*C)*(vd^2);
gsvd=2*m2*A*B*vd;


end

