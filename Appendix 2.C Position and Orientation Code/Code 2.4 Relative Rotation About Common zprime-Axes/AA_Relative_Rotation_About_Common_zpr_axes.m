%Given unit vectors f1pr,g1pr, and f2pr, Euler parameters p1 and p2, and 
%Euler parameter derivatives p1d and p2d for bodies 1 and 2 with common
%z1 and z2 axes, evaluate relative rotation and angular velocity

%Data are for Numerical Example
f1pr=[1;0;0;];
g1pr=[0;1;0];
f2pr=[1;0;0];
p1=[1;0;0;0];
e02=cos(pi/4);
e2=sin(pi/4)*[0;0;1];
p2=[e02;e2];

p1d=-0.5*GT(p1)'*[0;0;1];
p2d=0.5*GT(p2)'*[0;0;1];

%Compute theta12 in Eq. (3.4.37)
A1=AT(p1);
A2=AT(p2);
f1=A1*f1pr;
g1=A1*g1pr;
f2=A2*f2pr;

s=g1'*f2;
c=f1'*f2;
Arcs=asin(s);       %asin(s) is Arcsin(s)

if s<=0&&c<0
theta12=-pi-Arcs;
end 
if c>=0
theta12=Arcs;
end
if s>=0&&c<0
theta12=p1-Arcs;
end

%Compute theta12d in Eq. (3.4.40)
A1=AT(p1);
A2=AT(p2);

f1=A1*f1pr;
g1=A1*g1pr;
f2=A2*f2pr;

s=g1'*f2;
c=f1'*f2;

cgmsf1pr=c*g1pr-s*f1pr;
cgmsf1=c*g1-s*f1;

theta12d=f2'*BT(p1,cgmsf1pr)*p1d+cgmsf1'*BT(p2,f2pr)*p2d;
