function P3=P3(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

% Enter Third Derivative of Constraints P3=((((Psq)qd)sq)qd)sq

[r1,phi1,r2,phi2]=qPart(q);
[r1d,phi1d,r2d,phi2d]=qdPart(qd);

P3=[zeros(1,2),-(phi1d^2)*sin(phi1),zeros(1,3);
    zeros(1,2),(phi1d^2)*cos(phi1),zeros(1,3);
    zeros(1,2),(phi1d^2)*sin(phi1),zeros(1,2),(phi2d^2)*sin(phi2);
    zeros(1,2),-(phi1d^2)*cos(phi1),zeros(1,2),-(phi2d^2)*cos(phi2)];


end

