function QA=QAEval(t,q,qd,Lam,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

[q1,q2,q3]=qpart(q);
[q1d,q2d,q3d]=qpart(qd);

[Sfr1,Sfrpr1] = SfrSfrpr(q1d,par);
[Sfr2,Sfrpr2] = SfrSfrpr(q2d,par);
[Sfr3,Sfrpr3] = SfrSfrpr(q3d,par);

[cs2,dcs2]=csign(Lam*q2,par);
[cs1,dcs1]=csign((Lam*q1-m2*g),par);

QA=[-m1*g-K1*q1-Lam*q2*cs2*Sfr1;...
    K2*(q3-q2-1)-(Lam*q1-m2*g)*cs1*Sfr2;...
    -K2*(q3-q2-1)-m3*g*Sfr3];
end
    



