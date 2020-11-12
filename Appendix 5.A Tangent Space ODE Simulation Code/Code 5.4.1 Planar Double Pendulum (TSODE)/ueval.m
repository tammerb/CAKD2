function [u,Iteru] = ueval(t,ue,v,q0,V,U,B,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

err=utol+1;
i=1;
q=q0+V*v-U*ue;
u=ue;
P=Phi(t,q,par);
while err>utol    
delu=B*P;
u=u+delu;
q=q0+V*v-U*u;
P=Phi(t,q,par);
i=i+1;
err=norm(P);
end
Iteru=i-1;
end

