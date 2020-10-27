%Tangent Space RK DAE Exp; 3-Particle Friction-Piecewise Coulomb

utol=0.00000001;     %Tolerance in solving for u
Btol=0.00000001;     %Convergence criteria in B iteration
intol=0.000001;     %Tolerance in solving discretized  quations of motion
Maxv=0.7;       %Limit on magnitude of v
Maxiter=6;      %Limit on number of integration iterations

h=0.0001;
tfinal=3.99;

mode=1; %(Set to 1 to start)mode=1, no stic; mode=2, q1 const; 
                            %mode=3, q2 const; mode=4, q3 const
SticConst=2;   %SticConst=1 for q1 constraint,  2 for q2 constraint 
SticVelTol=0.001;   %Slip Velocity stiction tolerance
SticCritTol=0.001;    %Slip Criteria tolerance

g=9.8;
K1=10;   %Spring Constants
K2=10;
m1=6.5;
m2=2;
m3=2;
el=5;   %Length of distance constraint

mud=0.3;  %Coefficient of Dynamic Friction
mus=0.5;

nq=3;
nh1=1;
nv1=2;
nu1=1;
nh2=2;
nv2=1;
nu2=2;

q1b=0;
q2b=0;
q3b=0;

%Enter all parameters to be used
par=[nq;nh1;nv1;nu1;nh2;nv2;nu2;g;m1;m2;m3;K1;K2;el;...
    mud;mus;utol;Btol;intol;h];

Integ=1;  %Integ=1, Trap PartNR; %Integ=2, Trap FullNR;

% Data Storage Arrays
Qd=zeros(nq,10);
Q=zeros(nq,10);
Qdd=zeros(nq,10);
Vv1=zeros(nv1,10);
Vv1d=zeros(nv1,10);
Vv1dd=zeros(nv1,10);
Uu1=zeros(nu1,10);
Uu1d=zeros(nu1,10);
Uu1dd=zeros(nu1,10);
LLam1=zeros(nh1,10);
Vv2=zeros(nv2,10);
Vv2d=zeros(nv2,10);
Vv2dd=zeros(nv2,10);
Uu2=zeros(nu2,10);
Uu2d=zeros(nu2,10);
Uu2dd=zeros(nu2,10);
LLam2=zeros(nh2,10);

% Initial Conditions
t(1)=0;
q0=[0;5;6];
qd0=[0;0;-1];

if mode==1
%Initial qdd and Lam

    Lam0=0;

Phiq=PhiqEval(q0,par,mode);
Gam=GamEval(q0,qd0,par,mode);

QA=QAEval(q0,qd0,Lam0,par,mode);
AccCoef=[MEval(q0,par),Phiq';Phiq,0];
AccRHS=[QA;-Gam];

z=AccCoef\AccRHS;

AccCoefcond=cond(AccCoef);

qdd0=[z(1);z(2);z(3)];
Lam0=z(4);

Q(:,1)=q0;
Qd(:,1)=qd0;
Qdd(:,1)=qdd0;
LLam1(:,1)=Lam0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);
q=q0;
qd=qd0;
qdd=qdd0;

RMcond(1)=0;
end

%Start Integration Process
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
R1nrpt(1)=0;
SolIter(1)=Maxiter+1; 
jRepar=0;       %Counter for Reparameterization
Cr=2;
jiterrpt(1)=0;
vnormrpt(1)=Maxv+1;
Sw21=0;
Sw12=0;
npar=1;
q0=zeros(nq,1);
U1=zeros(nq,nu1);
V1=zeros(nq,nv1);
B1=zeros(nu1,nu1);
U2=zeros(nq,nu2);
V2=zeros(nq,nv2);
B2=zeros(nu2,nu2);
n0=1;
while t(n)<tfinal

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

jiter=jiterrpt(n-1);

if vnormrpt(n-1)>Maxv
Cr=2;
end

if jiter>Maxiter
Cr=Cr+2;
end
   
Crrpt(n)=Cr;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mode==1
  
[q,qd,qdd,Lam1,v1,v1d,v1dd,u1,q0,U1,V1,B1,jRepar,npar,...
    Rnorm,CondJ,Sw21,jiter,Cr]=...
INTEGRATE1(n,Q,Qd,Qdd,Vv1,Vv1d,Vv1dd,LLam1,Uu1,q0,U1,V1,B1,par,jRepar,...
npar,Cr,mode,q1b,q2b,q3b,Sw21);

EqErr(n)=Rnorm;
CondJ(n)=CondJ;

%Process Results
Vv1(:,n)=v1;
Vv1d(:,n)=v1d;
Vv1dd(:,n)=v1dd;
Uu1(:,n)=u1;
LLam1(:,n)=Lam1;
Vv2(:,n)=zeros(nv2,1);
Vv2d(:,n)=zeros(nv2,1);
Vv2dd(:,n)=zeros(nv2,1);
Uu2(:,n)=zeros(nu2,1);
LLam2(:,n)=zeros(nu2,1);
vnormrpt(n)=norm(v1);
vdnormrpt(n)=norm(v1d);
vddnormrpt(n)=norm(v1dd);
Lamnormrpt(n)=norm(Lam1);
Q(:,n)=q;
qnormrpt(n)=norm(q);
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

%Calculate Total Energy
M=MEval(q,par);
TE(n)=0.5*qd'*M*qd+m1*g*q(1)+0.5*K2*(q(3)-q(2)-1)^2+0.5*K1*q(1)^2;

Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

jiterrpt(n)=jiter;

%Friction Forces

fr1(n)=-mud*abs(Lam1*q(1))*sign(qd(1));
fr2(n)=-mud*abs(Lam1*q(1)-m2*g)*sign(qd(2));

end
%End mode ==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mode>1
    
[q,qd,qdd,Lam2,v2,v2d,v2dd,u2,q0,U2,V2,B2,jRepar,npar,Rnorm,...
    CondJ,Sw12,jiter,Cr]=...
    INTEGRATE2(n,Q,Qd,Qdd,Vv2,Vv2d,Vv2dd,LLam2,Uu2,q0,U2,V2,B2,...
    par,jRepar,npar,Cr,mode,Sw12,q1b,q2b,q3b);

EqErr(n)=Rnorm;
CondJ(n)=CondJ;

%Process Results
Vv2(:,n)=v2;
Vv2d(:,n)=v2d;
Vv2dd(:,n)=v2dd;
Uu2(:,n)=u2;
LLam2(:,n)=Lam2;
Vv1(:,n)=zeros(nv1,1);
Vv1d(:,n)=zeros(nv1,1);
Vv1dd(:,n)=zeros(nv1,1);
Uu1(:,n)=zeros(nu1,1);
LLam1(:,n)=zeros(nu1,1);
vnormrpt(n)=norm(v2);
vdnormrpt(n)=norm(v2d);
vddnormrpt(n)=norm(v2dd);
Lamnormrpt(n)=norm(Lam2);
Q(:,n)=q;
qnormrpt(n)=norm(q);
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

%Calculate Total Energy
M=MEval(q,par);
TE(n)=0.5*qd'*M*qd+m1*g*q(1)+0.5*K2*(q(3)-q(2)-1)^2+0.5*K1*q(1)^2;

Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);
jiterrpt(n)=jiter;


end
%End Mode > 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data of Interest
y(n)=q(1);
x1(n)=q(2);
x2(n)=q(3);
yd(n)=qd(1);
x1d(n)=qd(2);
x2d(n)=qd(3);
ydd(n)=qdd(1);
x1dd(n)=qdd(2);
x2dd(n)=qdd(3);

qd12=norm([qd(1);qd(2)]);
qd12norm(n)=qd12;
qdd12norm(n)=norm([qdd(1);qdd(2)]);

%Calculate constraint error
PosConstrNorm(n) = norm(PhiEval(q,par,mode,q1b,q2b,q3b));
Phiq=PhiqEval(q,par,mode);
VelConstrNorm(n)=norm(Phiq*qd);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par,mode);
AccelErr=Phiq*qdd+Gam;
AccConstrNorm(n)=norm(AccelErr);

%Stiction Criteria For Masses 1 and 2

alpha=mus*q(2)-q(1);
alpharpt(n)=alpha;
beta=mus*q(2)+q(1);
betarpt(n)=beta;
gama=m1*g+K1*q(1);
delta=mus*q(1)+q(2);
deltarpt(n)=delta;
epsilon=mus*q(1)-q(2);
epsilonrpt(n)=epsilon;
eta=K2*(q(3)-q(2)-1)+mus*m2*g;
tau=K2*(q(3)-q(2)-1)-mus*m2*g;
min1=min([gama/alpha;-gama/beta]);
max1=max([gama/alpha;-gama/beta]);
min2=min([eta/delta;-tau/epsilon]);
max2=max([eta/delta;-tau/epsilon]);

if abs(min1)>100
    min1=100*sign(min1);
end
if abs(min2)>100
    min2=100*sign(min2);
end
if abs(max1)>100
    max1=100*sign(max1);
end
if abs(max2)>100
    max2=100*sign(max2);
end

Stict=0;    %Stict=0 implies no stiction possible

if(alpha*beta<0)&&(delta*epsilon<0)
part(n)=1;  %part is designation of partitions 1-4

maxmin=max([min1;min2]);
minmax=min([max1,max2]);
if maxmin<=minmax+SticCritTol
Stict=1;    %Stict=1 implies stiction possible

if q(1)>=0
sf1min(n)=-minmax*q(1)-m1*g-K1*q(1);
sf1max(n)=-maxmin*q(1)-m1*g-K1*q(1);
end
if q(1)<0
sf1min(n)=-maxmin*q(1)-m1*g-K1*q(1);
sf1max(n)=-minmax*q(1)-m1*g-K1*q(1);
end

if q(2)>=0
sf2min(n)=-minmax*q(2)+K2*(q(3)-q(2)-1);
sf2max(n)=-maxmin*q(2)+K2*(q(3)-q(2)-1);
end
if q(2)<0
sf2min(n)=-maxmin*q(2)+K2*(q(3)-q(2)-1);
sf2max(n)=-minmax*q(2)+K2*(q(3)-q(2)-1);
end

end
end
%end of part 1

if(alpha*beta>0)&&(delta*epsilon<0)
part(n)=2;

if min2<=min1+SticCritTol
Stict=1;

if q(1)>=0
sf1min(n)=-min1*q(1)-m1*g-K1*q(1);
sf1max(n)=-min2*q(1)-m1*g-K1*q(1);
end
if q(1)<0
sf1min(n)=-min2*q(1)-m1*g-K1*q(1);
sf1max(n)=-min1*q(1)-m1*g-K1*q(1);
end

if q(2)>=0
sf2min(n)=-min1*q(2)+K2*(q(3)-q(2)-1);
sf2max(n)=-min2*q(2)+K2*(q(3)-q(2)-1);
end
if q(2)<0
sf2min(n)=-min2*q(2)+K2*(q(3)-q(2)-1);
sf2max(n)=-min1*q(2)+K2*(q(3)-q(2)-1);
end

end

if max2>=max1-SticCritTol
Stict=1;

if q(1)>=0
sf1min(n)=-max2*q(1)-m1*g-K1*q(1);
sf1max(n)=-max1*q(1)-m1*g-K1*q(1);
end
if q(1)<0
sf1min(n)=-max1*q(1)-m1*g-K1*q(1);
sf1max(n)=-max2*q(1)-m1*g-K1*q(1);
end

if q(2)>=0
sf2min(n)=-max2*q(2)+K2*(q(3)-q(2)-1);
sf2max(n)=-max1*q(2)+K2*(q(3)-q(2)-1);
end
if q(2)<0
sf2min(n)=-max1*q(2)+K2*(q(3)-q(2)-1);
sf2max(n)=-max2*q(2)+K2*(q(3)-q(2)-1);
end

end
end
%end of part 2

if(alpha*beta<0)&&(delta*epsilon>0)
part(n)=3;

if min2>=min1-SticCritTol
Stict=1;
if q(1)>=0
    sf1min(n)=-min2*q(1)-m1-K1*q(1);
    sf1max(n)=-min1*q(1)-m1-K1*q(1);
end
if q(1)<0
    sf1min(n)=-min1*q(1)-m1-K1*q(1);
    sf1max(n)=-min2*q(1)-m1-K1*q(1);
end

if q(2)>=0
    sf2min(n)=-min2*q(2)+K2*(q(3)-q(2)-1);
    sf2max(n)=-min1*q(2)+K2*(q(3)-q(2)-1);
end
if q(2)<0
    sf1min(n)=-min1*q(2)+K2*(q(3)-q(2)-1);
    sf1max(n)=-min2*q(2)+K2*(q(3)-q(2)-1);
end
end

if max2<=max1+SticCritTol
Stict=1;
if q(1)>=0
    sf1min(n)=-max1*q(1)-m1-K1*q(1);
    sf1max(n)=-max2*q(1)-m1-K1*q(1);
end
if q(1)<0
    sf1min(n)=-max2*q(1)-m1-K1*q(1);
    sf1max(n)=-max1*q(1)-m1-K1*q(1);
end

if q(2)>=0
    sf2min(n)=-max1*q(2)+K2*(q(3)-q(2)-1);
    sf2max(n)=-max2*q(2)+K2*(q(3)-q(2)-1);
end
if q(2)<0
    sf2min(n)=-max2*q(2)+K2*(q(3)-q(2)-1);
    sf2max(n)=-max1*q(2)+K2*(q(3)-q(2)-1);
end

end
end
%end part 3

if(alpha*beta>0)&&(delta*epsilon>0)
part(n)=4;
Stict=1;
end
%end part 4

Stiction(n)=Stict;

%mode change criteria and action

if mode==1
    
    %criteria for q1 or q2  stick
    SlpCrit12=max([abs(qd(1));abs(qd(2))]); %Slip criteria for qd1&qd2
    if (SlpCrit12<SticVelTol)&&(Stict==1)
        
        if SticConst==2
        q2b=q(2);   %value of stuck q2
        mode=2;
        Cr=2;
        Sw12=2; %Switching from 1 constraint to 2 constraints
        end
        
        if SticConst==1
        q1b=q(1);   %value of stuck q1
        mode=3;
        Cr=2;
        Sw12=2; %Switching from 1 constraint to 2 constraints
        end
    end
    %criteria for q3  stick
    q3frcoef=K2*(q(3)-q(2)-1)/(m2*g);
    SlpCrit3=abs(qd(3));
    if (SlpCrit3<SticVelTol)&&(abs(K2*(q(3)-q(2)-1))<=mus*m2*g) 
        mode=4;
        q3b=q(3);   %value of stuck q3
        Cr=2;
        Sw12=2; %Switching from 1 constraint to 2 constraints
    end
end

if mode==2
    if Sw12<1
    if Stict==0
        mode=1;
        Cr=2;
        Sw21=2; %Switching from 2 constraints to 1 constraint
    end
    end
end

if mode==3
    if Sw12<1
    if Stict==0
        mode=1;
        Cr=2;
        Sw21=2; %Switching from 2 constraints to 1 constraint
    end
    end
end

if mode==4
    if Sw12<1
        mode4FrictCoef=(abs(Lam2(2))/(m3*g));
        mode4FrictCoefrpt(n)=mode4FrictCoef;
    if mode4FrictCoef>mus
        mode=1;
        Cr=2;
        Sw21=2; %Switching from 2 constraints to 1 constraint
    end
    end
end

    Sw12rpt(n)=Sw12;
    moderpt(n)=mode;
    
    SlpCrit12rpt(n)=max([abs(qd(1));abs(qd(2))]);
    



end


   