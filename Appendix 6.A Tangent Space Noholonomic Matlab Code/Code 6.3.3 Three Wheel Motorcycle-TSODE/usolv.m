function [u,uiter]=usolv(t,ue,v,q0,V,U,B,utol,par,J,dpP0,dpP1,dpP2,...
    bp,ux,uy,uz,P,A1,apppsa,atpppsa)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
err=utol+1;
i=1;
q=q0+V*v-U*ue;
u=ue;
while err>utol;
delu=B*P0(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(delu);
end
uiter=i-1;
end

