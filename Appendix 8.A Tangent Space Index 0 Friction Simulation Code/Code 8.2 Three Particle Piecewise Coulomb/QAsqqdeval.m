function [QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

if mode==1
M=MEval(q,par);
QAsq=[-K1,-mud*sign(Lam*q(2))*sign(qd(1)),0;...
    -mud*Lam*sign(Lam*q(1)-m2*g)*sign(qd(2)),-K2,K2;0,K2,-K2];
QAsqd=zeros(3,3);
QAsLam=[-mud*q(2)*sign(Lam*q(2))*sign(qd(1));...
    -mud*q(1)*sign(Lam*q(1)-m2*g)*sign(qd(2));0];
end

if mode==2
M=MEval(q,par);
QAsq=[-K1,0,0;...
    0,-K2,K2;0,K2,-K2];
QAsqd=zeros(3,3);
QAsLam=zeros(3,2);
end

if mode==3
M=MEval(q,par);
QAsq=[-K1,0,0;...
    0,-K2,K2;0,K2,-K2];
QAsqd=zeros(3,3);
QAsLam=zeros(3,2);
end

if mode==4
M=MEval(q,par);
QAsq=[-K1,-mud*sign(Lam(1)*q(2))*sign(qd(1)),0;...
    -mud*Lam(1)*sign(Lam(1)*q(1)-m2*g)*sign(qd(2)),-K2,K2;0,K2,-K2];
QAsqd=zeros(3,3);
QAsLam=[-mud*q(2)*sign(Lam(1)*q(2))*sign(qd(1)),0;...
    -mud*q(1)*sign(Lam(1)*q(1)-m2*g)*sign(qd(2)),0;0,0];

end

