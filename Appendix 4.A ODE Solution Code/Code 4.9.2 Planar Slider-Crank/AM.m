function M=AM(v,par,dat)

% Enter Mass Matrix M=M(v,par)
[nv,intol,Atol,hmax,hvar]=BparPart(par);
[R,J1,m2,omega0]=AdatPart(dat);

s=sin(v);
c=cos(v);
A=(R*s+(R^2)*s*c/sqrt(4-(R*s)^2));

M=J1+m2*A^2;


end

