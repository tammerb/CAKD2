function Phiq=P1(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
%Enter Jacobian of holonomic constraint functions; example is tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[thet,thetd,thetdd]=Steer(t,par);
p2=[cos(thet/2);uz*sin(thet/2)];
AT=ATran(p);
A2=ATran(p2);
appp=[0;a];
atppp=[0;P*a];
c32=uz'*(BTran(p,bp)-s*BTran(p,A1*uz)+BTran(p,A1*A2*appp));
Phiq=[uz',uz'*BTran(p,dpP1),0,0,0;uz',uz'*BTran(p,dpP2),0,0,0;...
    uz',c32,uz'*AT*A1*A2*apppsa,-uz'*AT*A1*uz;...
    0,0,0,uz'*BTran(p,A1*A2*atppp),uz'*AT*A1*A2*atpppsa,0;...
    0,0,0,p',0,0,0;0,0,0,0,0,0,0,a',0];
end
