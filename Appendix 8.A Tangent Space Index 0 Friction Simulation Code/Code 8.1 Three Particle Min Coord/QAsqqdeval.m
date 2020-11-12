function [QAsq,QAsqd,QAsLam]=QAsqqdeval(t,q,qd,Lam,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

[q1,q2,q3]=qpart(q);
[q1d,q2d,q3d]=qpart(qd);

[Sfr1,Sfrpr1] = SfrSfrpr(q1d,par);
[Sfr2,Sfrpr2] = SfrSfrpr(q2d,par);
[Sfr3,Sfrpr3] = SfrSfrpr(q3d,par);

[cs2,dcs2]=csign(Lam*q2,par);
[cs1,dcs1]=csign((Lam*q1-m2*g),par);

QAsq=[-K1,-(Lam*cs2+(Lam^2)*q2*dcs2)*Sfr1,0;...
    -(Lam*cs1+Lam*(Lam*q1-m2*g)*dcs1)*Sfr2,-K2,K2;0,K2,-K2];
QAsqd=[-Lam*q2*cs2*Sfrpr1,0,0;0,-(Lam*q1-m2*g)*cs1*Sfrpr2,0;...
    0,0,-m3*g*Sfrpr3];
QAsLam=[-(q2*cs2+Lam*(q2^2)*dcs2)*Sfr1;...
    -(q1*cs1+q1*(Lam*q1-m2*g)*dcs1)*Sfr2;0];


    
end




