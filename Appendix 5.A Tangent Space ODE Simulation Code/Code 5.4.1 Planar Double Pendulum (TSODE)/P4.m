function P4=P4(t,q,n,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

% Enter P4=(PsqT)nsq

[r1,phi1,r2,phi2]=qPart(q);

P4=[zeros(2,nq);zeros(1,2),(n(1)-n(3))*cos(phi1)+(n(2)-n(4))*sin(phi1),...
    zeros(1,3);zeros(2,nq);zeros(1,5),-(n(3)+n(4))*sin(phi2)];


end

