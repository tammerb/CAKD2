function [u,Iteru] = ueval(t,ue,v,q0,V,U,B,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

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

