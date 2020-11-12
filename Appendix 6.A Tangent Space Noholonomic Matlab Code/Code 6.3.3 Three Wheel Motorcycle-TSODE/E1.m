function E=E1(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
%Enter coefficient matrix of differential constraint; example is tricycle
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[thet,thetd,thetdd]=Steer(t,par);
p2=[cos(thet/2);uz*sin(thet/2)];
AT=ATran(p);
A2=ATran(p2);
appp=[0;a];
atppp=[0;P*a];
at=AT*A1*A2*atppp;
atuz=atil(uz);
E=[ux'*AT',ux'*AT'*BTran(p,dpP0),0,0,0;...
    at'*atuz,at'*atuz*(BTran(p,bp)-s*BTran(p,A1*uz)+...
    BTran(p,A1*A2*appp)),at'*atuz*AT*A1*A2*apppsa,-at'*atuz*AT*A1*uz];

end

