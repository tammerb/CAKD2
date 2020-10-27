function P3=P3Eval(tn,q,qd,SJDT,par)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA]=parPart(par);

P3=zeros(nc,ngc);

I3=eye(3);
k=1;        %Joint No.
m=0;        %Constraint Counter-1
while k<=nh

%Distance Constraint
if SJDT(1,k)==1
[i,j,s1pr,s2pr,d]=DistPart(k,SJDT);

[P31,P32]=bbP3dist(i,j,s1pr,s2pr,d,tn,q,qd,par);
P3=Add(P3,P31,m,7*(i-1));      

if j>=1
P3=Add(P3,P32,m,7*(j-1));
end

m=m+1;
end

%Spherical Constraint
if SJDT(1,k)==2
[i,j,s1pr,s2pr]=SphPart(k,SJDT);
%P3 spherical contribution is zero, no nonzero terms to add

m=m+3;
end

%Cylindrical Constraint
if SJDT(1,k)==3
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=CylPart(k,SJDT);    
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

[P3cyl11,P3cyl12]=bbP3dot2(i,j,ux2pr,s1pr,s2pr,tn,q,qd,par);
[P3cyl21,P3cyl22]=bbP3dot2(i,j,uy2pr,s1pr,s2pr,tn,q,qd,par);
[P3cyl31,P3cyl32]=bbP3dot1(i,j,uz1pr,ux2pr,tn,q,qd,par);
[P3cyl41,P3cyl42]=bbP3dot1(i,j,uz1pr,uy2pr,tn,q,qd,par);

P31k=[P3cyl11;P3cyl21;P3cyl31;P3cyl41];
P32k=[P3cyl12;P3cyl22;P3cyl32;P3cyl42];

P3=Add(P3,P31k,m,7*(i-1));

if j>=1
P3=Add(P3,P32k,m,7*(j-1));
end

m=m+4;
end

%Revolute Constraint
if SJDT(1,k)==4
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=RevPart(k,SJDT);    
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

[P3cyl11,P3cyl12]=bbP3dot2(i,j,ux2pr,s1pr,s2pr,tn,q,qd,par);
[P3cyl21,P3cyl22]=bbP3dot2(i,j,uy2pr,s1pr,s2pr,tn,q,qd,par);
[P3cyl31,P3cyl32]=bbP3dot1(i,j,uz1pr,ux2pr,tn,q,qd,par);
[P3cyl41,P3cyl42]=bbP3dot1(i,j,uz1pr,uy2pr,tn,q,qd,par);
[P3rev51,P3rev52]=bbP3dot2(i,j,uz2pr,s1pr,s2pr,tn,q,qd,par);

P31k=[P3cyl11;P3cyl21;P3cyl31;P3cyl41;P3rev51];
P32k=[P3cyl12;P3cyl22;P3cyl32;P3cyl42;P3rev52];

P3=Add(P3,P31k,m,7*(i-1));

if j>=1
P3=Add(P3,P32k,m,7*(j-1));
end

m=m+5;
end

%Translational Constraint
if SJDT(1,k)==5
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=TranPart(k,SJDT);    
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

[P3cyl11,P3cyl12]=bbP3dot2(i,j,ux2pr,s1pr,s2pr,tn,q,qd,par);
[P3cyl21,P3cyl22]=bbP3dot2(i,j,uy2pr,s1pr,s2pr,tn,q,qd,par);
[P3cyl31,P3cyl32]=bbP3dot1(i,j,uz1pr,ux2pr,tn,q,qd,par);
[P3cyl41,P3cyl42]=bbP3dot1(i,j,uz1pr,uy2pr,tn,q,qd,par);
[P3tran51,P3tran52]=bbP3dot1(i,j,uy1pr,ux2pr,tn,q,qd,par);

P31k=[P3cyl11;P3cyl21;P3cyl31;P3cyl41;P3tran51];
P32k=[P3cyl12;P3cyl22;P3cyl32;P3cyl42;P3tran52];

P3=Add(P3,P31k,m,7*(i-1));

if j>=1
P3=Add(P3,P32k,m,7*(j-1));
end

m=m+5;
end

%Universal Constraint
if SJDT(1,k)==6
[i,j,s1pr,s2pr,uz1pr,uz2pr]=UnivPart(k,SJDT);    
[P3ks1,P3ks2]=bbP3sph(i,j,s1pr,s2pr,tn,q,par);
[P3kd11,P3kd12]=bbP3dot1(i,j,uz1pr,uz2pr,tn,q,par);
P3k1=[P3ks1;P3kd11];
P3k2=[P3ks2;P3kd12];
P3=Add(P3,P3k1,m,7*(i-1));

if j>=1
P3=Add(P3,P3k2,m,7*(j-1));
end
m=m+4;
end

%Strut Constraint
if SJDT(1,k)==7
[i,j,s1pr,s2pr,ux1pr,uy1pr]=StrutPart(k,SJDT);
[P3str11,P3str12]=bbP3dot2(i,j,ux1pr,s1pr,s2pr,tn,q,par);
[P3str21,P3str22]=bbP3dot2(i,j,uy1pr,s1pr,s2pr,tn,q,par);
P3k1=[P3str11;P3str21];
P3k2=[P3str12;P3str22];
P3=Add(P3,P3k1,m,7*(i-1));

if j>=1
P3=Add(P3,P3k2,m,7*(j-1));
end
m=m+2;
end

%Revolute-Spherical Constraint
if SJDT(1,k)==8
[i,j,s1pr,s2pr,d,uz2pr]=RevSphPart(k,SJDT);
[P3RS11,P3RS12]=bbP3dist(i,j,s1pr,s2pr,d,tn,q,par);
[P3RS21,P3RS22]=bbP3dot2(i,j,uz2pr,s1pr,s2pr,tn,q,par);
P3k1=[P3RS11;P3RS21];
P3k2=[P3RS12;P3RS22];
P3=Add(P3,P3k1,m,7*(i-1));

if j>=1
P3=Add(P3,P3k2,m,7*(j-1));
end
m=m+2;
end

k=k+1;

end

%Euler Parameter Normalization Constraint
%P3 Euler parameter contribution is zero, no nonzero terms to add


end





