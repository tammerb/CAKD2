function M=Meval(q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
%Enter Mass Matrix nqxnq; example below is for tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q)
I3=eye(3);
G=Geval(p);
M=blkdiag(m*I3,4*G'*J*G,zeros(3,3));


end

