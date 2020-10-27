function QA=QAEval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,...
    P,A1,apppsa,atpppsa)
%Evaluate Generalized Appoied Force nqx1; example below is for tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[rd,pd,ad,sd]=qdPart(qd)
AT=ATran(p);
G=Geval(p);
Fs=-k*(s-sf)-c*sd;
QA=[-m*g*uz+F*AT*uy;2*F*G'*ux;0;0;Fs];


end

