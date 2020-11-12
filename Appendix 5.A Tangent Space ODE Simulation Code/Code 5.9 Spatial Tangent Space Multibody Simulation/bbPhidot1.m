function Phi=bbPhidot1(i,j,a1pr,a2pr,tn,q,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
Phi=a1pr'*A1'*a2pr;    
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
Phi=a1pr'*A1'*A2*a2pr;
end   
        
end
