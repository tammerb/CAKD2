function P2=P2Eval(tn,q,x,SJDT,par)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA]=parPart(par);

P2=zeros(nc,ngc);

I3=eye(3);
k=1;        %Joint No.
m=0;        %Constraint Counter-1
while k<=nh

%Distance Constraint
if SJDT(1,k)==1
[i,j,s1pr,s2pr,d]=DistPart(k,SJDT);

[P21,P22]=bbP2dist(i,j,s1pr,s2pr,d,tn,q,x,par);
P2=Add(P2,P21,m,7*(i-1));

if j>=1
P2=Add(P2,P22,m,7*(j-1));
end

m=m+1;
end

%Spherical Constraint
if SJDT(1,k)==2
[i,j,s1pr,s2pr]=SphPart(k,SJDT);

[P21,P22]=bbP2sph(i,j,s1pr,s2pr,tn,q,x,par);
P2=Add(P2,P21,m,7*(i-1));        

if j>=1    
P2=Add(P2,P22,m,7*(j-1));
end

m=m+3;
end

%Cylindrical Constraint
if SJDT(1,k)==3
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=CylPart(k,SJDT);    
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

[P2cyl11,P2cyl12]=bbP2dot2(i,j,ux2pr,s1pr,s2pr,tn,q,x,par);
[P2cyl21,P2cyl22]=bbP2dot2(i,j,uy2pr,s1pr,s2pr,tn,q,x,par);
[P2cyl31,P2cyl32]=bbP2dot1(i,j,uz1pr,ux2pr,tn,q,x,par);
[P2cyl41,P2cyl42]=bbP2dot1(i,j,uz1pr,uy2pr,tn,q,x,par);

P21k=[P2cyl11;P2cyl21;P2cyl31;P2cyl41];
P22k=[P2cyl12;P2cyl22;P2cyl32;P2cyl42];

P2=Add(P2,P21k,m,7*(i-1));

if j>=1
P2=Add(P2,P22k,m,7*(j-1));
end

m=m+4;

end

%Revolute Constraint
if SJDT(1,k)==4
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=RevPart(k,SJDT);    
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

[P2cyl11,P2cyl12]=bbP2dot2(i,j,ux2pr,s1pr,s2pr,tn,q,x,par);
[P2cyl21,P2cyl22]=bbP2dot2(i,j,uy2pr,s1pr,s2pr,tn,q,x,par);
[P2cyl31,P2cyl32]=bbP2dot1(i,j,uz1pr,ux2pr,tn,q,x,par);
[P2cyl41,P2cyl42]=bbP2dot1(i,j,uz1pr,uy2pr,tn,q,x,par);
[P2rev51,P2rev52]=bbP2dot2(i,j,uz2pr,s1pr,s2pr,tn,q,x,par);

P21k=[P2cyl11;P2cyl21;P2cyl31;P2cyl41;P2rev51];
P22k=[P2cyl12;P2cyl22;P2cyl32;P2cyl42;P2rev52];

P2=Add(P2,P21k,m,7*(i-1));

if j>=1
P2=Add(P2,P22k,m,7*(j-1));
end

m=m+5;
end

%Translational Constraint
if SJDT(1,k)==5
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=TranPart(k,SJDT);    
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

[P2cyl11,P2cyl12]=bbP2dot2(i,j,ux2pr,s1pr,s2pr,tn,q,x,par);
[P2cyl21,P2cyl22]=bbP2dot2(i,j,uy2pr,s1pr,s2pr,tn,q,x,par);
[P2cyl31,P2cyl32]=bbP2dot1(i,j,uz1pr,ux2pr,tn,q,x,par);
[P2cyl41,P2cyl42]=bbP2dot1(i,j,uz1pr,uy2pr,tn,q,x,par);
[P2tran51,P2tran52]=bbP2dot1(i,j,uy1pr,ux2pr,tn,q,x,par);

P21k=[P2cyl11;P2cyl21;P2cyl31;P2cyl41;P2tran51];
P22k=[P2cyl12;P2cyl22;P2cyl32;P2cyl42;P2tran52];

P2=Add(P2,P21k,m,7*(i-1));

if j>=1
P2=Add(P2,P22k,m,7*(j-1));
end

m=m+5;
end

%Universal Constraint
if SJDT(1,k)==6
[i,j,s1pr,s2pr,uz1pr,uz2pr]=UnivPart(k,SJDT);    
[P2ks1,P2ks2]=bbP2sph(i,j,s1pr,s2pr,tn,q,par);
[P2kd11,P2kd12]=bbP2dot1(i,j,uz1pr,uz2pr,tn,q,par);
P2k1=[P2ks1;P2kd11];
P2k2=[P2ks2;P2kd12];
P2=Add(P2,P2k1,m,7*(i-1));

if j>=1
P2=Add(P2,P2k2,m,7*(j-1));
end
m=m+4;
end

%Strut Constraint
if SJDT(1,k)==7
[i,j,s1pr,s2pr,ux1pr,uy1pr]=StrutPart(k,SJDT);
[P2str11,P2str12]=bbP2dot2(i,j,ux1pr,s1pr,s2pr,tn,q,par);
[P2str21,P2str22]=bbP2dot2(i,j,uy1pr,s1pr,s2pr,tn,q,par);
P2k1=[P2str11;P2str21];
P2k2=[P2str12;P2str22];
P2=Add(P2,P2k1,m,7*(i-1));

if j>=1
P2=Add(P2,P2k2,m,7*(j-1));
end
m=m+2;
end

%Revolute-Spherical Constraint
if SJDT(1,k)==8
[i,j,s1pr,s2pr,d,uz2pr]=RevSphPart(k,SJDT);
[P2RS11,P2RS12]=bbP2dist(i,j,s1pr,s2pr,d,tn,q,par);
[P2RS21,P2RS22]=bbP2dot2(i,j,uz2pr,s1pr,s2pr,tn,q,par);
P2k1=[P2RS11;P2RS21];
P2k2=[P2RS12;P2RS22];
P2=Add(P2,P2k1,m,7*(i-1));

if j>=1
P2=Add(P2,P2k2,m,7*(j-1));
end
m=m+2;
end

k=k+1;

end

%Euler Parameter Normalization Constraints
i=1;
while i<=nb
[xr1,xp1]=xPart(x,i);
P21=[zeros(1,3),xp1'];
P2=Add(P2,P21,m,7*(i-1));
i=i+1;
m=m+1;
end

end





