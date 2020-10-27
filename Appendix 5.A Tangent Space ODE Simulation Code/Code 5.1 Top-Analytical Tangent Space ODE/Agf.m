function gev=Agf(t,v,vd,par,p0,V,Jpr,m,g,Dcf)

%Enter right side of ODE
p=p0+V*v-(1-sqrt(1-v'*v))*p0;
D=(eye(4)-p0*p'/(p'*p0))*V;
Gp=Geval(p);
pd=D*vd;
Gpd=Geval(pd);
Ap=BAeval(p);
uz=[0;0;1];

gev=4*(pd'*pd/(p'*p0))*D'*Gp'*Jpr*Gp*p0+8*D'*Gpd'*Jpr*Gpd*p-...
    2*D'*Gp'*(m*g*Batil(uz)*Ap'*uz+4*Dcf*pd'*Gp'*Gp*pd*uz);

end

