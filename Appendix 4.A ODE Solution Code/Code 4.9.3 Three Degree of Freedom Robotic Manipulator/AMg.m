function [M,gf,M2,gfsv,gfsvd] = AMg(t,v,vd,vdd,par,dat,integ,deriv)

[nv,intol,Atol,hmax,hvar]=BparPart(par);
[g,m1,m2,m3,J1ii,J21,J22,J23,K1,K2,K3]=AdatPart(dat);

m1=200;
m2=100;
m3=100;
J1=J1ii*eye(3);
J2=diag([J21;J22;J23]);
v1=v(1);
v2=v(2);
v3=v(3);
v1d=vd(1);
v2d=vd(2);
v3d=vd(3);
c1=cos(v1);
s1=sin(v1);
c2=cos(v2);
s2=sin(v2);
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

A1=[c1,-s1,0;s1,c1,0;0,0,1];
A12=[1,0,0;0,c2,-s2;0,s2,c2];
A1p=[-s1,-c1,0;c1,-s1,0;0,0,0];
A12p=[0,0,0;0,-s2,-c2;0,c2,-s2];
A1pp=[-c1,s1,0;-s1,-c1,0;0,0,0];
A12pp=[0,0,0;0,-c2,s2;0,-s2,-c2];

W1=[zeros(2,3);1,0,0];
V2=[A1p*A12*uy,A1*A12p*uy,zeros(3,1)];
W2=[A12'*A1'*uz,ux,zeros(3,1)];
V3=[(v3+2)*A1p*A12*uy,(v3+2)*A1*A12p*uy,A1*A12*uy];

U2=A1pp*A12*uy*v1d^2+2*A1p*A12p*uy*v1d*v2d+A1*A12pp*uy*v2d^2;
X2=A12'*A1p'*uz*v1d^2+A12p'*A1'*uz*v1d*v2d;
U3=(v3+2)*A1pp*A12*uy*v1d^2+(v3+2)*A1*A12pp*uy*v2d^2+...
    2*(v3+2)*A1p*A12p*uy*v1d*v2d+2*A1p*A12*uy*v1d*v3d+2*A1*A12p*uy*v2d*v3d;

%Evaluate mass matrix
M=W1'*J1*W1+m2*V2'*V2+m3*V3'*V3+2*W2'*J2*W2;

%Evaluate right side of ODE
QA=0*sin(t)*[1;1;1];
QI=[-K1*v1;-K2*v2;-K3*v3];

gf=-(W1'*atil(W1*vd)*J1*W1*vd+m2*V2'*(U2+g*uz)+m3*V3'*(U3+g*uz)...
    +2*W2'*J2*X2+2*W2'*atil(W2*vd)*J2*W2*vd)+QA+QI;

if integ<5
%Set M2,gfsv,and gfsvd to zero for explicit integration
Z=zeros(3,3);
M2=Z;
gfsv=Z;
gfsvd=Z;
end

if integ>4
%Implicit integration
if deriv==1
%Retain M and gf for residual calculation and zeros for derivatives    
Z=zeros(3,3);
M2=Z;
gfsv=Z;
gfsvd=Z;
end

if deriv==2
%Calculate M2,gfsv,and gfsvd  

a=-A1p*A12*uy*v1d^2+2*A1pp*A12p*uy*v1d*v2d+A1p*A12pp*uy*v2d^2;
b=A1pp*A12p*uy*v1d^2+2*A1p*A12pp*uy*v1d*v2d-A1*A12p*uy&v2d^2;
U2sv=[a,b,zeros(3,1)];

U2svd=[2*A1pp*A12*uy*v1d+2*A1p*A12p*uy*v2d,2*A1p*A12p*uy*v1d+...
    2*A1*A12pp*uy*v2d,zeros(3,1)];

c=-(v3+2)*A1p*A12*uy*v1d^2+(v3+2)*A1pp*A12p*uy*v1d*v2d+...
    (v3+2)*A1p*A12pp*uy*v2d^2+2*A1pp*A12*uy*v1d*v3d+2*A1p*A12p*uy*v2d*v3d;
d=(v3+2)*A1pp*A12p*uy*v1d^2+(v3+2)*A1p*A12pp*uy*v1d*v2d-...
    (v3+2)*A1*A12p*uy*v2d^2+2*A1p*A12p*uy*v1d*v3d+2*A1*A12pp*uy*v2d*v3d;    
e=A1pp*A12*uy*v1d^2+2*A1p*A12p*uy*v1d*v2d+A1*A12pp*uy*v2d^2;
U3sv=[c,d,e];

f2=2*(v3+2)*A1pp*A12*uy*v1d+2*(v3+2)*A1p*A12p*uy*v2d+2*A1p*A12*uy*v3d;
g2=2*(v3+2)*A1p*A12p*uy*v1d+2*(v3+2)*A1*A12pp*uy*v2d+2*A1*A12p*uy*v3d;
h2=2*A1p*A12*uy*v1d+2*A1*A12p*uy*v2d;
U3svd=[f2,g2,h2];

alph=vdd;
V2alphsv=[alph(1)*A1pp*A12*uy+alph(2)*A1p*A12p*uy,...
    alph(1)*A1p*A12p*uy+alph(2)*A1*A12pp*uy,zeros(3,1)];

alph=vdd;
i=alph(1)*(v3+2)*A1pp*A12*uy+alph(2)*(v3+2)*A1p*A12p*uy+alph(3)*A1p*A12*uy;
j=alph(1)*(v3+2)*A1p*A12p*uy+alph(2)*(v3+2)*A1*A12pp*uy+alph(3)*A1*A12p*uy;
k=alph(1)*A1p*A12*uy+alph(2)*A1*A12p*uy;
V3alphsv=[i,j,k];

bet=U2+g*uz;
V2Tbetsv=[uy'*A12'*A1pp'*bet,uy'*A12p'*A1p'*bet,0;...
    uy'*A12p'*A1p'*bet,uy'*A12pp'*A1'*bet,0;0,0,0];

bet=U3+g*uz;
V3Tbetsv=[(v3+2)*uy'*A12'*A1pp'*bet,(v3+2)*uy'*A12p'*A1p'*bet,...
    uy'*A12'*A1p'*bet;...
    (v3+2)*uy'*A12p'*A1p'*bet,(v3+2)*uy'*A12pp'*A1'*bet,...
    uy'*A12p'*A1'*bet;uy'*A12'*A1p'*bet,uy'*A12p'*A1'*bet,0];

X2sv=[A12'*A1pp'*uz*v1d^2+A12p'*A1p'*uz*v1d*v2d,...
    A12p'*A1p'*uz*v1d^2+A12pp'*A1'*uz*v1d*v2d,zeros(3,1)];

X2svd=[2*A12'*A1p'*uz*v1d+A12p'*A1'*uz*v2d,A12p'*A1'*uz*v2d,zeros(3,1)];

gam=vdd;
W2gamsv1=[gam(1)*A12'*A1p'*uz,gam(1)*A12p'*A1'*uz,zeros(3,1)];
gam=vd;
W2gamsv2=[gam(1)*A12'*A1p'*uz,gam(1)*A12p'*A1'*uz,zeros(3,1)];

del=atil(W2*vd)*J2*W2*vd;
W2Tdelsv1=[uz'*A1p*A12*del,uz'*A1*A12p*del,0;zeros(2,3)];
del=J2*X2;
W2Tdelsv2=[uz'*A1p*A12*del,uz'*A1*A12p*del,0;zeros(2,3)];


M2=2*(2*m2*V2'*V2alphsv+2*W2'*J2*W2gamsv1+2*m3*V3'*V3Tbetsv);

gfsv=-(m2*V2Tbetsv+m2*V2'*U2sv+W2Tdelsv2+W2'*J2*X2sv+2*W2Tdelsv1+...
    2*W2'*(atil(W2*vd)*J2-atil(J2*W2*vd))*W2gamsv2+...
    m3*V3Tbetsv+m3*V3'*U3sv);

gfsvd=-(W1'*(atil(W1*vd)*J1-atil(J1*W1*vd))*W1+m2*V2'*U2svd+...
    2*W2'*J2*X2svd+2*W2'*(atil(W2*vd)*J2-atil(J2*W2*vd))*W2+m3*V3'*U3svd);


end    %End Calculation of M2,gfsv,and gfsvd

end

