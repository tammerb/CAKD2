function [P,Pst,Pstt]=P5Eval(tn,q,SJDT,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

% Enter Constraint t derivatives of P(t); Default is Zeros

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app ==1  %2 Bar with 2 rotational drivers

%Rot driver body1 to ground
k=3;
[i,j,vx1pr,vy1pr,vx2pr]=RotDrPart(k,SJDT);

a1=1;
a=a1;
theta=-a*sin(tn);  %Angle from ground, bod 0, to bod 1 is negative of angle 
                   %from bod 1 to bod 0
thetad=-a*cos(tn);
thetadd=a*sin(tn);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

PD1=PD;
PDd1=PDd;
PDdd1=PDdd;

%Rot driver body1 to body 2
k=4;
[i,j,vx1pr,vy1pr,vx2pr]=RotDrPart(k,SJDT);

a2=1;
a=a2;
theta=a*sin(tn);
thetad=a*cos(tn);
thetadd=-a*sin(tn);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

PD2=PD;
PDd2=PDd;
PDdd2=PDdd;

P=[zeros(nhc,1);PD1;PD2;zeros(nb,1)];
Pst=[zeros(nhc,1);PDd1;PDd2;zeros(nb,1)];
Pstt=[zeros(nhc,1);PDdd1;PDdd2;zeros(nb,1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app ==2  %4 Bar with internal rotational driver
k=4;
[i,j,vx1pr,vy1pr,vx2pr]=RotDrPart(k,SJDT);
    
a=1;
theta=a*sin(tn);
thetad=a*cos(tn);
thetadd=-a*sin(tn);

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

P=[zeros(nhc,1);PD;zeros(nb,1)];
Pst=[zeros(nhc,1);PDd;zeros(nb,1)];
Pstt=[zeros(nhc,1);PDdd;zeros(nb,1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app ==3  %2-body Slider, x-z plane, distDr Bod1-grnd
a=1;
b=3;
%d=b+a*sin(tn)
PD=-0.5*(b+a*sin(tn))^2;
P=[zeros(nhc,1);PD;zeros(nb,1)];
PDd=-(b+a*sin(tn))*a*cos(tn);
Pst=[zeros(nhc,1);PDd;zeros(nb,1)];
PDdd=(b+a*sin(tn))*a*sin(tn)-(a*cos(tn))^2;
Pstt=[zeros(nhc,1);PDdd;zeros(nb,1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app ==4  %Slider-Crank in x-z Plane, Bod 1 RotDr about x axis
k=4;
[i,j,vx1pr,vy1pr,vx2pr]=RotDrPart(k,SJDT);

omega=1;
theta=omega*tn;
thetad=omega;
thetadd=0;

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

P=[zeros(nhc,1);PD;zeros(nb,1)];
Pst=[zeros(nhc,1);PDd;zeros(nb,1)];
Pstt=[zeros(nhc,1);PDdd;zeros(nb,1)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app ==5  %Spatial Slider-Crank, RotDr Bod1-grnd
k=4;
[i,j,vx1pr,vy1pr,vx2pr]=RotDrPart(k,SJDT);

omega=4;
theta=omega*tn;
thetad=omega;
thetadd=0;

[r1,p1]=qPart(q,i);
A1=ATran(p1);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
PD=-sin(theta);
PDd=-thetad*cos(theta);
PDdd=-thetadd*cos(theta)+(thetad^2)*sin(theta);
else
PD=-cos(theta); 
PDd=thetad*sin(theta);
PDdd=thetadd*sin(theta)+(thetad^2)*cos(theta);
end   
end

P=[zeros(nhc,1);PD;zeros(nb,1)];
Pst=[zeros(nhc,1);PDd;zeros(nb,1)];
Pstt=[zeros(nhc,1);PDdd;zeros(nb,1)];
end

end