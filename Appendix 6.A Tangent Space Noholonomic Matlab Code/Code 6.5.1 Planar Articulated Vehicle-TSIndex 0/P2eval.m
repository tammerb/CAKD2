function P2=P2eval(t,q,x,par)

[r1,phi1,r2,phi2]=qPart(q);
[xr1,xphi1,xr2,xphi2]=xPart(x);
A1=ATran(phi1);
A2=ATran(phi2);
ipr=[1;0];

P2=[zeros(2,2),3*xphi1*A1*ipr,zeros(2,2),3*xphi2*A2*ipr];
end

