function M=MEval(q,SMDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

M=zeros(ngc,ngc);
I3=eye(3);
i=1;
while i<=nb
m=SMDT(1,i);
J=diag([SMDT(2,i);SMDT(3,i);SMDT(4,i)]);
[r,p]=qPart(q,i);
G=GEval(p);
Mi=blkdiag(m*I3,4*G'*J*G);
M=Add(M,Mi,7*(i-1),7*(i-1));
i=i+1;
end

end

