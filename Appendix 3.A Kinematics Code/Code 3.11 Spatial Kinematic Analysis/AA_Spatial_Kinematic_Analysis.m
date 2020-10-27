%AA Spatial Kinematic Analysis

qtol=10^-6;     %Tolerance in solving for u
h=0.01;       %Time step

tfinal=4*pi;

%Application Data
    %app=1, 2-Bar with 2 rotational drivers
    %app=2, 4-Bar with internal rotational driver
    %app=3, 2-body Slider in x-z plane, distDr Bod1-grnd
    %app=4, Slider-Crank in y-z Plane, Bod 1 RotDr about x axis
    %app=5, Spatial Slider-Crank, RotDr Bod1-grnd
    
app=5;

[nb,ngc,nh,nhc,nc,nd,SJDT,q0e]=AppData(app);

par=[nb;ngc;nh;nhc;nd;qtol;app];

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
z3=zeros(3,1);

% Data Storage Arrays
Q=zeros(ngc,10);
Qd=zeros(ngc,10);
Qdd=zeros(ngc,10);

Q(:,1)=q0e;

%Kinematic Analysis

n=1;
t(1)=0;
while t(n)<tfinal
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

%Evaluate q

q=Q(:,n-1);     %q-estimate
if n-1>1
q=q+h*Qd(:,n-1);
end

i=1;    %Position iteration counter
err=qtol+1;
while err >qtol
Phi=PhiEval(tn,q,SJDT,par);
Phiq=PhiqEval(q,SJDT,par);
if i==1
CondPhiq(n)=cond(Phiq);    
end
delq=-Phiq\Phi;
q=q+delq;
err=norm(Phi);
errrpt(i)=err;
i=i+1;
end
    
iter(n)=i-1;
Q(:,n)=q;

%Evaluate qd
Phiq=PhiqEval(q,SJDT,par);
[P,Pst,Pstt]=P5Eval(tn,q,SJDT,par);
qd=-Phiq\Pst;
Qd(:,n)=qd;

% Evaluate qdd
%P2=P2Eval(q,qd,SJDT,par);
%Gam=P2*qd+Pstt;
Gam=GamEval(tn,q,qd,SJDT,par);
qdd=-Phiq\Gam;
Qdd(:,n)=qdd;

%Calculate output data of interest (Enter for each application)

if app==1  %2 Bar with 2 rotational drivers 
thetaz1(1)=0;
thetaz2(1)=0;
thetaz2m1(1)=0;
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
[r2,p2]=qPart(q,2);
[r2d,p2d]=qPart(qd,2);
[r2dd,p2dd]=qPart(qdd,2);
omegaz1(n)=2*uz'*EEval(p1)*p1d;
thetaz1(n)=thetaz1(n-1)+h*omegaz1(n);
omegaz2(n)=2*uz'*EEval(p2)*p2d;
thetaz2(n)=thetaz2(n-1)+h*omegaz2(n);
omegaz2m1(n)=omegaz2(n)-omegaz1(n);
thetaz2m1(n)=thetaz2m1(n-1)+h*omegaz2m1(n);
omegazd1(n)=2*uz'*EEval(p1)*p1dd;
omegazd2(n)=2*uz'*EEval(p2)*p2dd;
omegazd2m1(n)=omegazd2(n)-omegazd1(n);
y2(n)=q(9);
y2d(n)=qd(9);
y2dd(n)=qdd(9);
x2(n)=q(8);
x2d(n)=qd(8);
x2dd(n)=qdd(8);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2       %4 Bar with internal rotational driver
thetaz2m1(1)=0;
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
[r2,p2]=qPart(q,2);
[r2d,p2d]=qPart(qd,2);
[r2dd,p2dd]=qPart(qdd,2);
omegaz2m1(n)=2*uz'*EEval(p2)*p2d-2*uz'*EEval(p1)*p1d;
thetaz2m1(n)=thetaz2m1(n-1)+h*omegaz2m1(n);
omegazd2m1(n)=2*uz'*EEval(p2)*p2dd-2*ux'*EEval(p1)*p1dd;
y2(n)=q(9);
y2d(n)=qd(9);
y2dd(n)=qdd(9);
x2(n)=q(8);
x2d(n)=qd(8);
x2dd(n)=qdd(8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %2-body Slider in x-z plane, distDr Bod1-grnd
x1(n)=q(1);
z2(n)=q(10);  
z2d(n)=qd(10);
z2dd(n)=qdd(10);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Slider-Crank in x-z Plane, Bod 1 RotDr about x axis
thetax(1)=0;
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
p1ddnorm(n)=norm(p1dd);
omega1x(n)=-2*ux'*EEval(p1)*p1d;  
thetax(n)=thetax(n-1)+h*omega1x(n);
omega1xd(n)=-2*ux'*EEval(p1)*p1dd;
y2(n)=q(9);
y2d(n)=qd(9);
y2dd(n)=qdd(9);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5       %Spatial Slider-Crank, RotDr Bod1-grnd
thetax(1)=0;
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
p1ddnorm(n)=norm(p1dd);
omega1x(n)=-2*ux'*EEval(p1)*p1d;     
thetax(n)=thetax(n-1)+h*omega1x(n);
omega1xd(n)=-2*ux'*EEval(p1)*p1dd;
x2(n)=q(8);
x2d(n)=qd(8);
x2dd(n)=qdd(8);

end

end





