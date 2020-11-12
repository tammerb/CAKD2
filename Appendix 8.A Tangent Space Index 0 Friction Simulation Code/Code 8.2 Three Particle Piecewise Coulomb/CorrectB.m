function [B,Biter]=CorrectB(q,B,U,par,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);
i=1;
Cnorm=Btol+1;

if mode==1
I=eye(nu1);
end
if mode>1
I=eye(nu2);
end

Phiq=PhiqEval(q,par,mode);
while Cnorm>Btol
i=i+1;
delB=-B*Phiq*U*B+B;
B=B+delB;
C=B*Phiq*U-I;
Cnorm=norm(C);
Biter=i-1;
end
end

