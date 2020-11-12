function Gam=Gameval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,...
    uz,P,A1,apppsa,atpppsa)
%Enter expressions for eqqd, et, Ptq, Ptt, and Gam; example is tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[rd,pd,ad,sd]=qdPart(qd);
[thet,thetd,thetdd]=Steer(t,par);
p2=[cos(thet/2);uz*sin(thet/2)];
p2d=(thetd/2)*[-sin(thet/2);uz*cos(thet/2)];
p2dd=(thetdd/2)*[-sin(thet/2);uz*cos(thet/2)]-...
    (thetd/2)^2*[cos(thet/2);uz*sin(thet/2)];
atuz=atil(uz);
AT=ATran(p);
A2=ATran(p2);
appp=[0;a];
atppp=[0;P*a];
at=AT*A1*A2*atppp;
atst=AT*A1*BTran(p2,atppp)*p2d;
atstt=AT*A1*(BTran(p2,atppp)*p2dd+BTran(p2d,atppp)*p2d);
Ptt=[0;0;uz'*AT*A1*(BTran(p2,appp)*p2dd+BTran(p2d,appp)*p2d);...
    uz'*atstt;0;0];
Ptqqd=[0;0;uz'*(BTran(pd,A1*BTran(p2,appp)*p2d)*p+...
    AT*A1*MTran(p2,p2d)*apppsa*ad);...
    uz'*(BTran(pd,A1*BTran(p2,atppp)*p2d)*p+...
    AT*A1*MTran(p2,p2d)*atpppsa*ad);0;0];

Etqd=[0;atst'*atuz*(rd+(BTran(p,bp)-s*BTran(p,A1*uz)+...
    BTran(p,A1*A2*appp))*pd)+atst'*atuz*AT*A1*A2*apppsa*ad+...
    at'*atuz*MTran(p,pd)*A1*BTran(p2,appp)*p2d-...
    sd*atst'*atuz*AT*A1*uz+at'*atuz*AT*A1*BTran(p2,apppsa*ad)*p2d];
et=[0;-atst'*atil(uz)*AT*A1*BTran(p2,appp)*p2d-...
    at'*atuz*AT*A1*(BTran(p2,appp)*p2dd+BTran(p2d,appp)*p2d)];
eqqd=[0;-atst'*atil(uz)*(BTran(p,A1*BTran(p2,appp)*p2d)*pd+...
    AT*A1*MTran(p2,p2d)*apppsa*ad)+(AT*A1*BTran(p2,appp)*p2d)'*...
    atuz*(BTran(p,A1*A2*atppp)*pd+AT*A1*A2*atpppsa*ad)];
P2=P2eval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,...
    P,A1,apppsa,atpppsa);
E2=E2eval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
Gam=[P2+2*Ptqqd+Ptt;E2+Etqd-eqqd-et];



end



