function P3=P3eval(t,q,qd,par)

[r1,phi1,r2,phi2]=qPart(q);
[r1d,phi1d,r2d,phi2d]=qdPart(qd);
P=[0,-1;1,0];
A1=ATran(phi1);
A2=ATran(phi2);
ipr=[1;0];

P3=[zeros(2,2),3*(phi1d^2)*P*A1*ipr,zeros(2,2),3*(phi2d^2)*P*A2*ipr];

end

