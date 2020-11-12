function Phi=bbPhiRotDr(i,j,vx1pr,vy1pr,vx2pr,q,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
Phi=s;   
else
Phi=c;    
end   
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
Phi=s;   
else
Phi=c;    
end
end   
        
end
