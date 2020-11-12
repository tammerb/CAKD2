function Phi = Phi(t,q,par)

[p1,p2,r2] = qpart(q);
uz=[0;0;1];
r21=r2+Aeval(p2)*uz+Aeval(p1)*uz;
Phi=0.5*[p1'*p1-1;p2'*p2-1;r21'*r21-1];
end

