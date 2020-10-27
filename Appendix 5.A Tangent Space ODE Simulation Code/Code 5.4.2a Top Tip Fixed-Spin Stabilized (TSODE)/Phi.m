function Phi = Phi(t,q,par)

[r,p]=qPart(q);
uz=[0;0;1];

Phi=[r-Aeval(p)*uz;(p'*p-1)/2];
end

