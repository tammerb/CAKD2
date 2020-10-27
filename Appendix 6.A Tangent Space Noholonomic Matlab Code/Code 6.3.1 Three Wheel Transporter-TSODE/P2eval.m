function P2=P2eval(t,q,qd,par)
%Evaluate (phiq*x)q*qd, planar art veh
[r1,phi1,r2,phi2]=qPart(q);
[r1d,phi1d,r2d,phi2d]=qdPart(qd)
[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ]=Partpar(par);
P=[0,-1;1,0];
A1=ATran(phi1);
A2=ATran(phi2);
ipr=[1;0];

P2=10*(phi1^2)*A1*ipr+10*(phi2^2)*A2*ipr
end

