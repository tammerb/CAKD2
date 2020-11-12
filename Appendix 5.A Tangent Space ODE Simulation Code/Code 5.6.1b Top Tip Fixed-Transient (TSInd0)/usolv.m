function [u,Iteru] = usolv(u,v,q0,V,U,B,utol)
err=utol+1;
i=1;
q=q0+V*v-U*u;
while err>utol;
delu=B*PhiEval(q);
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(PhiEval(q));
end
Iteru=i-1;
end

