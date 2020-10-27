function Phi=bbPhidist(i,j,s1pr,s2pr,d,q,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
d12=s2pr-r1-A1*s1pr;
Phi=(d12'*d12-d^2)/2;    
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
d12=r2+A2*s2pr-r1-A1*s1pr;
Phi=(d12'*d12-d^2)/2;
end   
        
end


    





