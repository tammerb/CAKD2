function [u,Iteru]=usolv(tn,u,v,q0,SJDT,V,U,B,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

err=utol+1;
i=1;
q=q0+V*v-U*u;
while err>utol
Phi=PhiEval(tn,q,SJDT,par);
delu=B*Phi;
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(Phi);
end
Iteru=i-1;
end

