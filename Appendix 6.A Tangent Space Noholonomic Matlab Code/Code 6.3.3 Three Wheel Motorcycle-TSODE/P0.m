function Phi=P0(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
%Evaluate Constraint Function, nhx1; example is tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[thet,thetd,thetdd]=Steer(t,par);
p2=[cos(thet/2);uz*sin(thet/2)];
AT=ATran(p);
A2=ATran(p2);
appp=[0;a];
atppp=[0;P*a];
Phi=[uz'*(r+AT*dpP1)-0.3;uz'*(r+AT*dpP2)-0.3;...
    uz'*(r+AT*bp-s*AT*A1*uz+AT*A1*A2*appp);uz'*AT*A1*A2*atppp;...
    (p'*p-1)/2;(a'*a-0.09)/2];
end

