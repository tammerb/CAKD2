function Phi = Phi(t,q,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);

uz=[0;0;1];
AT=ATran(p);
Phi=[uz'*(r-AT*uz);(p'*p-1)/2];
end

