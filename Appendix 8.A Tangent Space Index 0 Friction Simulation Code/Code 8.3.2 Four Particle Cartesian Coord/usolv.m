function [u,Iteru] = usolv(u,v,q0,V,U,B,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

err=utol+1;
i=1;
q=q0+V*v-U*u;
while err>utol
delu=B*PhiEval(q,par);
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(PhiEval(q,par));
end
Iteru=i-1;
end

