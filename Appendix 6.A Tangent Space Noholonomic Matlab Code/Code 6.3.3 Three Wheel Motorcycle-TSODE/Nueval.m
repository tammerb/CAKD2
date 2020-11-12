function Nu=Nueval(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,...
    apppsa,atpppsa)
%jEnter expressions for e, Pt, and Nu; example is tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[thet,thetd,thetdd]=Steer(t,par);
p2=[cos(thet/2);uz*sin(thet/2)];
p2d=(thetd/2)*[-sin(thet/2);uz*cos(thet/2)];
AT=ATran(p);
A2=ATran(p2);
appp=[0;a];
atppp=[0;P*a];
at=AT*A1*A2*atppp;
e=[0;-at'*atil(uz)*AT*A1*BTran(p2,appp)*p2d];
Pt=[0;0;uz'*AT*A1*BTran(p2,appp)*p2d;uz'*AT*A1*BTran(p2,atppp)*p2d;0;0];
Nu=[-Pt;e];


end

