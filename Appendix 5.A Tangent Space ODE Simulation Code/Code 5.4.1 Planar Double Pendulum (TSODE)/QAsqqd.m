function [QAsq,QAsqd]=QAsqqd(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

[r1,phi1,r2,phi2]=qPart(q);
x1=r1(1);

% Evaluate Derivatives QAsq(t,q,qd,par) and QAsqd(t,q,qd,par)
QAsq=[zeros(2,6);zeros(1,2),-(K1+K2),zeros(1,2),K2;...
    zeros(2,6);zeros(1,2),K2,zeros(1,2),-K2];
QAsqd=[zeros(2,6);zeros(1,2),-(C1+C2),zeros(1,2),C2;...
    zeros(2,6);zeros(1,2),C2,zeros(1,2),-C2];   

end

