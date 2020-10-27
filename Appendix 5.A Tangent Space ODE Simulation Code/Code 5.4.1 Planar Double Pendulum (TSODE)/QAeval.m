function QA=QAeval(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

[r1,phi1,r2,phi2]=qPart(q);
x1=r1(1);
[r1d,phi1d,r2d,phi2d]=qdPart(qd);

QA=[0;-g;-(-n1+K1*(phi1+pi/2)+C1*phi1d-...
    K2*(phi2-phi1)-C2*(phi2d-phi1d));0;-g;...
    -(-n2+K2*(phi2-phi1)+C2*(phi2d-phi1d))];
end

