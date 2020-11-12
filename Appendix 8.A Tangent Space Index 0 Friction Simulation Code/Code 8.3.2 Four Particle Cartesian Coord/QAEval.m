function QA=QAEval(q,qd,Lam,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

[csa1,dcsa1]=csign(Lam(1),par);
[csa2,dcsa2]=csign(Lam(2),par);
[csa3,dcsa3]=csign(Lam(3),par);
[csa4,dcsa4]=csign(Lam(4),par);
[csa5,dcsa5]=csign(Lam(5),par);
[csa6,dcsa6]=csign(Lam(6),par);

[Sfr1,Sfrpr1]=SfrSfrpr(qd(1),par);
[Sfr5,Sfrpr5]=SfrSfrpr(qd(5),par);
[Sfr9,Sfrpr9]=SfrSfrpr(qd(9),par);
[Sfr10,Sfrpr10]=SfrSfrpr(qd(10),par);

if FR==1 % Trans Joint    
Fn1=Lam(1)*csa1+Lam(2)*csa2;
Fn5=Lam(3)*csa3+Lam(4)*csa4;
Fn9=Lam(5)*csa5+Lam(6)*csa6;
Fn10=m4*g;
end

if FR==2 % Cyl Joint    
Fn1=((Lam(1)^2+Lam(2)^2)^0.5);
Fn5=((Lam(3)^2+Lam(4)^2)^0.5);
Fn9=((Lam(5)^2+Lam(6)^2)^0.5);
Fn10=m4*g;
end


QA=[K4*(q(10)-q(1)-1)-K1*q(1)-Fn1*Sfr1;0;-m1*g;0;-K2*q(5)-Fn5*Sfr5;...
    -m2*g;0;0;...
    -m3*g-K3*q(9)-Fn9*Sfr9;-K4*(q(10)-q(1)-1)-Fn10*Sfr10];

end


