function QA=QAEval(q,qd,Lam,par,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);


if mode==1
QA=[-m1*g-K1*q(1)-mud*abs(Lam*q(2))*sign(qd(1));...
    K2*(q(3)-q(2)-1)-mud*abs(Lam*q(1)-m2*g)*sign(qd(2));...
    -K2*(q(3)-q(2)-1)-mud*m3*g*sign(qd(3))];
end

if mode==2
QA=[-m1*g-K1*q(1);...
    K2*(q(3)-q(2)-1);...
    -K2*(q(3)-q(2)-1)-mud*m3*g*sign(qd(3))];
end

if mode==3
QA=[-m1*g-K1*q(1);...
    K2*(q(3)-q(2)-1);...
    -K2*(q(3)-q(2)-1)-mud*m3*g*sign(qd(3))];
end

if mode==4
QA=[-m1*g-K1*q(1)-mud*abs(Lam(1)*q(2))*sign(qd(1));...
    K2*(q(3)-q(2)-1)-mud*abs(Lam(1)*q(1)-m2*g)*sign(qd(2));...
    -K2*(q(3)-q(2)-1)];
end

end

