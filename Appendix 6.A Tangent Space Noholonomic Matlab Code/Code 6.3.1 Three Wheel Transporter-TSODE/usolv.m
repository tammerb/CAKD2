function [u,uiter] = usolv(t,ue,v,q0,V,U,B,utol,par)
err=utol+1;
i=1;
q=q0+V*v-U*ue;
u=ue;
while err>utol;
delu=B*P0(t,q,par);
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(delu);
end
uiter=i-1;
end

