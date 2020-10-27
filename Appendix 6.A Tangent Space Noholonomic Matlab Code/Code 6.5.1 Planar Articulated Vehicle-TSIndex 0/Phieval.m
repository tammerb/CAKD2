function Phi=Phieval(t,q,par)

[r1,phi1,r2,phi2]=qPart(q);
A1=ATran(phi1);
A2=ATran(phi2);
ipr=[1;0];

Phi=[r1-3*A1*ipr-r2-3*A2*ipr];



end