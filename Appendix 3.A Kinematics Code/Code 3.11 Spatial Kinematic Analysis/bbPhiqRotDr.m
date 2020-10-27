function [Phiq1,Phiq2]=bbPhiqRotDr(i,j,vx1pr,vy1pr,vx2pr,q,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
Phiq1=[0,0,0,vx2pr'*BTran(p1,vy1pr)];   
else
Phiq1=[0,0,0,vx2pr'*BTran(p1,vx1pr)];    
end
Phiq2=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
Phiq1=[0,0,0,vx2pr'*A2'*BTran(p1,vy1pr)];
Phiq2=[0,0,0,vy1pr'*A1'*BTran(p2,vx2pr)]; 
else
Phiq1=[0,0,0,vx2pr'*A2'*BTran(p1,vx1pr)];
Phiq2=[0,0,0,vx1pr'*A1'*BTran(p2,vx2pr)];    
end
end   
        
end

