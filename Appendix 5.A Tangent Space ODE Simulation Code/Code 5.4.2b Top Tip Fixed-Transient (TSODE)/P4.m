function P4=P4(t,q,n,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

% Enter P4=(PsqT)nsq
uz=[0;0;1];
nk=[n(1);n(2);n(3)];
np=n(4);

P4=[zeros(3,7);zeros(4,3),-Keval(uz,nk)+np*eye(4)];


end

