function [B,Biter]=CorrectB(q,B,U,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);
i=1;
Cnorm=Btol+1;
I=eye(nu);
Phiq=PhiqEval(q,par);
while Cnorm>Btol
i=i+1;
delB=-B*Phiq*U*B+B;
B=B+delB;
C=B*Phiq*U-I;
Cnorm=norm(C);
Biter=i-1;
end
end

