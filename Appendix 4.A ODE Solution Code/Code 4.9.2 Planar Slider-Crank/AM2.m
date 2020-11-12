function M2=AM2(v,mu,par)


% Enter M2=(M(v,par)mu)sv, mu same dimension as v
[nv,intol,Atol,hmax,hvar]=BparPart(par);
[R,J1,m2,omega0]=AdatPart(dat);

s=sin(v);
c=cos(v);
A=(R*s+(R^2)*s*c/sqrt(4-(R*s)^2));
B=(R*c+(R^2)*(c^2-s^2)/sqrt(4-(R*s)^2)+(R^4)*((s*c)^2)/(4-(R*s)^2)^1.5);
M2=mu*2*m2*A*B;


end

