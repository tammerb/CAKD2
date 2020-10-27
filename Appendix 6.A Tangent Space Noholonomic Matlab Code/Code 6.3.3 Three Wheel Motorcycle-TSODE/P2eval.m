function P2=P2eval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,...
    P,A1,apppsa,atpppsa)
%Evaluate (phiq*x)q*qd, example is tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[rd,pd,ad,sd]=qdPart(qd)
[thet,thetd,thetdd]=Steer(t,par);
p2=[cos(thet/2);uz*sin(thet/2)];
AT=ATran(p);
A2=ATran(p2);
appp=[0;a];
atppp=[0;P*a];
c3=uz'*(BTran(pd,bp)*pd-2*sd*BTran(pd,A1*uz)*p-s*BTran(pd,A1*uz)*pd+...
    BTran(pd,A1*A2*appp)*pd+MTran(p,pd)*A1*A2*apppsa*ad+...
    BTran(pd,A1*A2*apppsa*ad)*p);
c4=uz'*(BTran(pd,A1*A2*atppp)*pd+MTran(p,pd)*A1*A2*atpppsa*ad+...
    BTran(pd,A1*A2*atpppsa*ad)*p);
P2=[uz'*BTran(pd,dpP1)*pd;uz'*BTran(pd,dpP2)*pd;c3;c4;pd'*pd;ad'*ad];
end

