function Phiq=Phiqeval(t,q,par)
%Enter Jacobian of holonomic constraint functions; example is planar art
%veh
[r1,phi1,r2,phi2]=qPart(q);
P=[0,-1;1,0];
A1=ATran(phi1);
A2=ATran(phi2);
ipr=[1;0];

Phiq=[eye(2),-3*P*A1*ipr,-eye(2),-3*P*A2*ipr];
end
