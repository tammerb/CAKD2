function Phi=P0(t,q,par)
%Evaluate Constraint Function, nhx1; example below is for planar
%articulated vehicle
[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ]=Partpar(par);
[r1,phi1,r2,phi2]=qPart(q)[r,phi]=qPart(q);A1=ATran(phi1);
A2=ATran(phi2);
ipr=[1;0];
Phi=[r1-10*A1*ipr-r2-10*A2*ipr];



end