function Phiq = P1(t,q,par)

%Evaluate Constraint Jacobian
[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);

uz=[0;0;1];
BT=BTran(p,uz);
Phiq=[uz',-uz'*BT;0,0,0,p'];
end
