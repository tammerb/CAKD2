function [QAsq,QAsqd]=QAsqqd(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

% Evaluate Derivatives QAsq(t,q,qd,par) and QAsqd(t,q,qd,par)

QAsq=zeros(nq,nq);
QAsqd=zeros(nq,nq);

end

