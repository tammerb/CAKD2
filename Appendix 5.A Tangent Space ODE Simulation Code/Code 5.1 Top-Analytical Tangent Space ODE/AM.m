function M=AM(v,par,p0,V,Jpr,m,g)

% Enter Mass Matrix M=M(v,par,Jpr,m,g)
p=p0+V*v-(1-sqrt(1-v'*v))*p0;
D=(eye(4)-(1/(p'*p0))*p0*p')*V;
G=Geval(p);

M=4*D'*G'*Jpr*G*D;

end

