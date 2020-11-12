function QA=QAEval(q,qd,Lam,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

[csa1,dcsa1]=csign(Lam(1),par);
[csa2,dcsa2]=csign(Lam(2),par);
[Sfr2,Sfrpr2]=SfrSfrpr(qd(2),par);
[Sfr3,Sfrpr3]=SfrSfrpr(qd(3),par);
[Sfr5,Sfrpr5]=SfrSfrpr(qd(5),par);

QA=[0;-m1*g-K1*q(2)-Lam(1)*csa1*Sfr2;K2*(q(5)-q(3)-1)-Lam(2)*csa2*Sfr3;...
    -m2*g;-K2*(q(5)-q(3)-1)-m3*g*Sfr5];
end



