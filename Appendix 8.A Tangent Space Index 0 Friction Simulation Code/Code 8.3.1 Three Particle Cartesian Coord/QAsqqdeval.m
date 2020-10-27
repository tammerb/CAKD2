function [QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

[csa1,dcsa1]=csign(Lam(1),par);
[csa2,dcsa2]=csign(Lam(2),par);
[Sfr2,Sfrpr2]=SfrSfrpr(qd(2),par);
[Sfr3,Sfrpr3]=SfrSfrpr(qd(3),par);
[Sfr5,Sfrpr5]=SfrSfrpr(qd(5),par);

QAsq=[zeros(1,5);0,-K1,0,0,0;0,0,-K2,0,K2;zeros(1,5);0,0,K2,0,-K2];

QAsqd=[zeros(1,5);0,-Lam(1)*csa1*Sfrpr2,0,0,0;0,0,-Lam(2)*csa2*Sfrpr3,0,0;...
    zeros(1,5);0,0,0,0,-m3*g*Sfrpr5];

QAsLam=[zeros(1,3);-(csa1+Lam(1)*dcsa1)*Sfr2,0,0;...
    0,-(csa2+Lam(2)*dcsa2)*Sfr3,0;zeros(2,3)];
end




