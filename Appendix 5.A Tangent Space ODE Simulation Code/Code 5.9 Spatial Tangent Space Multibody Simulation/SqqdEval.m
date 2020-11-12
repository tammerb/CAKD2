function [Ssq,Ssqd]=SqqdEval(q,qd,SMDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

Ssq=zeros(ngc,ngc);
Ssqd=zeros(ngc,ngc);
I3=eye(3);
i=1;
while i<=nb
m=SMDT(1,i);
J=diag([SMDT(2,i);SMDT(3,i);SMDT(4,i)]);
[r,p]=qPart(q,i);
[rd,pd]=qPart(qd,i);
G=GEval(p);
Gd=GEval(pd);
Ssqi=[zeros(3,7);zeros(4,3),8*Gd'*J*Gd];
Ssq=Add(Ssq,Ssqi,7*(i-1),7*(i-1));
Ssqdi=[zeros(3,7);zeros(4,3),8*TEval(J*Gd*p)-8*Gd'*J*G];
Ssqd=Add(Ssqd,Ssqdi,7*(i-1),7*(i-1));
i=i+1;
end


end

