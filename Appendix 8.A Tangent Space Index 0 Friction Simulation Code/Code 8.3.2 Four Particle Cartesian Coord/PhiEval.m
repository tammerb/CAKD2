function Phi = PhiEval(q,par)
%Evaluate constraint function

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

Phi=[q(2);q(3);q(4);q(6);q(7);q(8);...
    ((q(6)-q(3))^2+(q(5)-q(2))^2+(q(4)-q(1))^2-el1^2)/2;
    ((q(9)-q(6))^2+(q(8)-q(5))^2+(q(7)-q(4))^2-el2^2)/2];
    end





