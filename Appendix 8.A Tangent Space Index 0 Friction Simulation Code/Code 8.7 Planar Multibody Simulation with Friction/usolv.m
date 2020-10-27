function [u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

err=utol+1;
i=1;
q=q0+V*v-U*u;
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
while err>utol
Phi=PhiEval(tn,q,PJDT,par);

delu=B*Phi;
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(Phi);
end
Iteru=i-1;
end

