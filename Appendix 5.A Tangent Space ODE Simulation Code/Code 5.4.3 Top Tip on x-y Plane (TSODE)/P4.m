function P4=P4(t,q,n,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);

% Enter P4=(PsqT)nsq
uz=[0;0;1];

P4=[zeros(3,nq);zeros(4,3),-n(1)*Keval(uz,uz)+n(2)*eye(4)];


end

