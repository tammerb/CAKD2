function P1 = P1(t,q,par)
%Evaluate Constraint Jacobian
[r,p]=qPart(q);
uz=[0;0;1];

P1=[eye(3),-BTran(p,uz);0,0,0,p'];
end
