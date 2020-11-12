function [QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par)

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

Fn1sLam1=(csa1+Lam(1)*dcsa1);
Fn1sLam2=(csa2+Lam(2)*dcsa2);
Fn5sLam3=(csa3+Lam(3)*dcsa3);
Fn5sLam4=(csa4+Lam(4)*dcsa4);
Fn9sLam5=(csa5+Lam(5)*dcsa5);
Fn9sLam6=(csa6+Lam(6)*dcsa6);

end

if FR==2 % Cyl Joint    
Fn1=((Lam(1)^2+Lam(2)^2)^0.5);
Fn5=((Lam(3)^2+Lam(4)^2)^0.5);
Fn9=((Lam(5)^2+Lam(6)^2)^0.5);
Fn10=m4*g;

Fn1sLam1=((Lam(1)^2+Lam(2)^2)^-0.5)*Lam(1);
Fn1sLam2=((Lam(1)^2+Lam(2)^2)^-0.5)*Lam(2);
Fn5sLam3=((Lam(3)^2+Lam(4)^2)^-0.5)*Lam(3);
Fn5sLam4=((Lam(3)^2+Lam(4)^2)^-0.5)*Lam(4);
Fn9sLam5=((Lam(5)^2+Lam(6)^2)^-0.5)*Lam(5);
Fn9sLam6=((Lam(5)^2+Lam(6)^2)^-0.5)*Lam(6);
end


QAsq=[-K4-K1,zeros(1,8),K4;zeros(3,10);zeros(1,4),-K2,zeros(1,5);...
    zeros(3,10);zeros(1,8),-K3,0;K4,zeros(1,8),-K4];

QAsqd=[-Fn1*Sfrpr1,zeros(1,9);zeros(3,10);...
    zeros(1,4),-Fn5*Sfrpr5,zeros(1,5);zeros(3,10);...
    zeros(1,8),-Fn9*Sfrpr9,0;zeros(1,9),-Fn10*Sfrpr10];

QAsLam=-[Fn1sLam1*Sfr1,Fn1sLam2*Sfr1,zeros(1,6);zeros(3,8);...
    0,0,Fn5sLam3*Sfr5,Fn5sLam4*Sfr5,zeros(1,4);zeros(3,8);...
    zeros(1,4),Fn9sLam5*Sfr9,Fn9sLam6*Sfr9,0,0;zeros(1,8)];
end




