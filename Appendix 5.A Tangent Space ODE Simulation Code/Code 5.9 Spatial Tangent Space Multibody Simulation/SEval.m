function S=SEval(q,qd,SMDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

S=zeros(ngc,1);
I3=eye(3);
i=1;
while i<=nb
m=SMDT(1,i);
J=diag([SMDT(2,i);SMDT(3,i);SMDT(4,i)]);
[r,p]=qPart(q,i);
[rd,pd]=qPart(qd,i);
G=GEval(p);
Gd=GEval(pd);
Si=[0;0;0;8*Gd'*J*Gd*p];
S=Add(S,Si,7*(i-1),0);
i=i+1;
end


end

