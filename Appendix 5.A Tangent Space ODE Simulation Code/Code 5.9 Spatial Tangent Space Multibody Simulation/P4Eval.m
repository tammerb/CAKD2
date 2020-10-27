function P4=P4Eval(tn,q,eta,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

P4=zeros(ngc,ngc);

I3=eye(3);
z3=zeros(3,1);

k=1;    %Joint number
m=0;    %Address in vector eta

while k<=nh

%Distance Constraint
if SJDT(1,k)==1
[i,j,s1pr,s2pr,d]=DistPart(k,SJDT);
etak=eta(m+1);

[P411,P412,P422]=bbP4dist(i,j,s1pr,s2pr,d,tn,q,etak,par);
P4=Add(P4,P411,7*(i-1),7*(i-1));

if j>=1
P421=P412';
P4=Add(P4,P412,7*(i-1),7*(j-1));
P4=Add(P4,P421,7*(j-1),7*(i-1));
P4=Add(P4,P422,7*(j-1),7*(j-1));
end

m=m+1;
end

%Spherical Constraint
if SJDT(1,k)==2
[i,j,s1pr,s2pr]=SphPart(k,SJDT);
etak=[eta(m+1);eta(m+2);eta(m+3)];

[P411,P412,P422]=bbP4sph(i,j,s1pr,s2pr,tn,q,etak,par);    
P4=Add(P4,P411,7*(i-1),7*(i-1));    

if j>=1
%P412 and P421 are zero
P4=Add(P4,P422,7*(j-1),7*(j-1));
end

m=m+3;
end

%Cylindrical Constraint
if SJDT(1,k)==3
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=CylPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

[P4cyl111,P4cyl112,P4cyl122]=...
    bbP4dot2(i,j,vx2pr,s1pr,s2pr,tn,q,eta(m+1),par);
[P4cyl211,P4cyl212,P4cyl222]=...
    bbP4dot2(i,j,vy2pr,s1pr,s2pr,tn,q,eta(m+2),par);
[P4cyl311,P4cyl312,P4cyl322]=bbP4dot1(i,j,vz1pr,vx2pr,tn,q,eta(m+3),par);
[P4cyl411,P4cyl412,P4cyl422]=bbP4dot1(i,j,vz1pr,vy2pr,tn,q,eta(m+4),par);

P411k=P4cyl111+P4cyl211+P4cyl311+P4cyl411;
P412k=P4cyl112+P4cyl212+P4cyl312+P4cyl412;
P421k=P412k';
P422k=P4cyl122+P4cyl222+P4cyl322+P4cyl422;

P4=Add(P4,P411k,7*(i-1),7*(i-1));

if j>=1
P4=Add(P4,P412k,7*(i-1),7*(j-1));
P4=Add(P4,P421k,7*(j-1),7*(i-1));
P4=Add(P4,P422k,7*(j-1),7*(j-1));
end

m=m+4;
end

%Revolute Constraint
if SJDT(1,k)==4
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=RevPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

[P4cyl111,P4cyl112,P4cyl122]=...
    bbP4dot2(i,j,vx2pr,s1pr,s2pr,tn,q,eta(m+1),par);
[P4cyl211,P4cyl212,P4cyl222]=...
    bbP4dot2(i,j,vy2pr,s1pr,s2pr,tn,q,eta(m+2),par);
[P4cyl311,P4cyl312,P4cyl322]=bbP4dot1(i,j,vz1pr,vx2pr,tn,q,eta(m+3),par);
[P4cyl411,P4cyl412,P4cyl422]=bbP4dot1(i,j,vz1pr,vy2pr,tn,q,eta(m+4),par);
[P4rev511,P4rev512,P4rev522]=...
    bbP4dot2(i,j,vz2pr,s1pr,s2pr,tn,q,eta(m+5),par);

P411k=P4cyl111+P4cyl211+P4cyl311+P4cyl411+P4rev511;
P412k=P4cyl112+P4cyl212+P4cyl312+P4cyl412+P4rev512;
P421k=P412k';
P422k=P4cyl122+P4cyl222+P4cyl322+P4cyl422+P4rev522;

P4=Add(P4,P411k,7*(i-1),7*(i-1));

if j>=1
P4=Add(P4,P412k,7*(i-1),7*(j-1));
P4=Add(P4,P421k,7*(j-1),7*(i-1));
P4=Add(P4,P422k,7*(j-1),7*(j-1));
end

m=m+5;
end

%Translational Constraint
if SJDT(1,k)==5
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=TranPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

[P4cyl111,P4cyl112,P4cyl122]=...
    bbP4dot2(i,j,vx2pr,s1pr,s2pr,tn,q,eta(m+1),par);
[P4cyl211,P4cyl212,P4cyl222]=...
    bbP4dot2(i,j,vy2pr,s1pr,s2pr,tn,q,eta(m+2),par);
[P4cyl311,P4cyl312,P4cyl322]=bbP4dot1(i,j,vz1pr,vx2pr,tn,q,eta(m+3),par);
[P4cyl411,P4cyl412,P4cyl422]=bbP4dot1(i,j,vz1pr,vy2pr,tn,q,eta(m+4),par);
[P4tran511,P4tran512,P4tran522]=...
    bbP4dot1(i,j,vy1pr,vx2pr,tn,q,eta(m+5),par);

P411k=P4cyl111+P4cyl211+P4cyl311+P4cyl411+P4tran511;
P412k=P4cyl112+P4cyl212+P4cyl312+P4cyl412+P4tran512;
P421k=P412k';
P422k=P4cyl122+P4cyl222+P4cyl322+P4cyl422+P4tran522;

P4=Add(P4,P411k,7*(i-1),7*(i-1));

if j>=1
P4=Add(P4,P412k,7*(i-1),7*(j-1));
P4=Add(P4,P421k,7*(j-1),7*(i-1));
P4=Add(P4,P422k,7*(j-1),7*(j-1));
end

m=m+5;
end

%Universal Constraint
if SJDT(1,k)==6
[i,j,s1pr,s2pr,vz1pr,vz2pr]=UnivPart(k,SJDT);    
[P4ks11,P4ks12,P4ks22]=bbP4sph(i,j,s1pr,s2pr,tn,q,par);
[P4kd111,P4kd112,P4kd122]=bbP4dot1(i,j,vz1pr,vz2pr,tn,q,par);
P4k11=P4ks11+P4kd111;
P4k12=P4ks12+P4kd112;
P4k22=P4ks22+P4kd122;
P4k21=P4k12';

P4=Add(P4,P4k11,7*(i-1),7*(i-1));

if j>=1
P4=Add(P4,P4k12,7*(i-1),7*(j-1)); 
P4=Add(P4,P4k21,7*(j-1),7*(i-1));
P4=Add(P4,P4k22,7*(j-1),7*(j-1));
end

m=m+4;
end

%Strut Constraint
if SJDT(1,k)==7
[i,j,s1pr,s2pr,vx1pr,vy1pr]=StrutPart(k,SJDT);
[P4str111,P4str112,P4str122]=bbP4dot2(i,j,vx1pr,s1pr,s2pr,tn,q,par);
[P4str211,P4str212,P4str222]=bbP4dot2(i,j,vy1pr,s1pr,s2pr,tn,q,par);
P4k11=P4str111+P4str211;
P4k12=P4str112+P4str212;
P4k22=P4str122+P4str222;
P4k21=P4k12';

P4=Add(P4,P4k11,7*(i-1),7*(i-1));

if j>=1
P4=Add(P4,P4k12,7*(i-1),7*(j-1)); 
P4=Add(P4,P4k21,7*(j-1),7*(i-1));
P4=Add(P4,P4k22,7*(j-1),7*(j-1));
end

m=m+2;
end

%Revolute-Spherical Constraint
if SJDT(1,k)==8
[i,j,s1pr,s2pr,d,vz2pr]=RevSphPart(k,SJDT);
[P4RS111,P4RS112,P4RS122]=bbP4dist(i,j,s1pr,s2pr,d,tn,q,par);
[P4RS211,P4RS212,P4RS222]=bbP4dot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);
P4k11=P4RS111+P4RS211;
P4k12=P4RS112+P4RS212;
P4k22=P4RS122+P4RS222;
P4k21=P4k12';

P4=Add(P4,P4k11,7*(i-1),7*(i-1));

if j>=1
P4=Add(P4,P4k12,7*(i-1),7*(j-1)); 
P4=Add(P4,P4k21,7*(j-1),7*(i-1));
P4=Add(P4,P4k22,7*(j-1),7*(j-1));
end

m=m+2;
end

k=k+1;
end

%Euler Parameter Normalization Constraint
i=1;
while i<=nb
etak=eta(m+1);
P411=etak*[zeros(3,7);zeros(4,3),eye(4)];
P4=Add(P4,P411,7*(i-1),7*(i-1));
i=i+1;
m=m+1;
end

end







