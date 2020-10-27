function Phi=PhiEval(tn,q,SJDT,par)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA]=parPart(par);

Phi=zeros(nc,1);
I3=eye(3);
k=1;        %Joint No.
m=0;        %Constraint Counter-1
while k<=nh

%Distance Constraint
if SJDT(1,k)==1
[i,j,s1pr,s2pr,d]=DistPart(k,SJDT);

Phik=bbPhidist(i,j,s1pr,s2pr,d,tn,q,par);
Phi=Add(Phi,Phik,m,0);
m=m+1;
end

%Spherical Constraint
if SJDT(1,k)==2
[i,j,s1pr,s2pr]=SphPart(k,SJDT);    
Phik=bbPhisph(i,j,s1pr,s2pr,tn,q,par);
Phi=Add(Phi,Phik,m,0);
m=m+3;
end

%Cylindrical Constraint
if SJDT(1,k)==3
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=CylPart(k,SJDT);
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

Phicyl1=bbPhidot2(i,j,ux2pr,s1pr,s2pr,tn,q,par);
Phicyl2=bbPhidot2(i,j,uy2pr,s1pr,s2pr,tn,q,par);
Phicyl3=bbPhidot1(i,j,uz1pr,ux2pr,tn,q,par);
Phicyl4=bbPhidot1(i,j,uz1pr,uy2pr,tn,q,par);

Phik=[Phicyl1;Phicyl2;Phicyl3;Phicyl4];
Phi=Add(Phi,Phik,m,0);
m=m+4;
end

%Revolute Constraint
if SJDT(1,k)==4
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=RevPart(k,SJDT);
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

Phicyl1=bbPhidot2(i,j,ux2pr,s1pr,s2pr,tn,q,par);
Phicyl2=bbPhidot2(i,j,uy2pr,s1pr,s2pr,tn,q,par);
Phicyl3=bbPhidot1(i,j,uz1pr,ux2pr,tn,q,par);
Phicyl4=bbPhidot1(i,j,uz1pr,uy2pr,tn,q,par);
Phirev5=bbPhidot2(i,j,uz2pr,s1pr,s2pr,tn,q,par);

Phik=[Phicyl1;Phicyl2;Phicyl3;Phicyl4;Phirev5];
Phi=Add(Phi,Phik,m,0);
m=m+5;
end

%Translational Constraint
if SJDT(1,k)==5
[i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=TranPart(k,SJDT);
uy1pr=atil(uz1pr)*ux1pr;   
uy2pr=atil(uz2pr)*ux2pr;

Phicyl1=bbPhidot2(i,j,ux2pr,s1pr,s2pr,tn,q,par);
Phicyl2=bbPhidot2(i,j,uy2pr,s1pr,s2pr,tn,q,par);
Phicyl3=bbPhidot1(i,j,uz1pr,ux2pr,tn,q,par);
Phicyl4=bbPhidot1(i,j,uz1pr,uy2pr,tn,q,par);
Phitran5=bbPhidot1(i,j,uy1pr,ux2pr,tn,q,par);

Phik=[Phicyl1;Phicyl2;Phicyl3;Phicyl4;Phitran5];
Phi=Add(Phi,Phik,m,0);
m=m+5;
end

%Universal Constraint
if SJDT(1,k)==6
[i,j,s1pr,s2pr,uz1pr,uz2pr]=UnivPart(k,SJDT);    
Phiks=bbPhisph(i,j,s1pr,s2pr,tn,q,par);
Phikd1=bbPhidot1(i,j,uz1pr,uz2pr,tn,q,par);
Phik=[Phiks;Phikd1];
Phi=Add(Phi,Phik,m,0);
m=m+4;
end

%Strut Constraint
if SJDT(1,k)==7
[i,j,s1pr,s2pr,ux2pr,uy2pr]=StrutPart(k,SJDT);
Phistr1=bbPhidot2(i,j,ux2pr,s1pr,s2pr,tn,q,par);
Phistr2=bbPhidot2(i,j,uy2pr,s1pr,s2pr,tn,q,par);
Phik=[Phistr1;Phistr2];
Phi=Add(Phi,Phik,m,0);
m=m+2;
end

%Revolute-Spherical Constraint
if SJDT(1,k)==8
[i,j,s1pr,s2pr,d,uz2pr]=RevSphPart(k,SJDT);
PhiRS1=bbPhidist(i,j,s1pr,s2pr,d,tn,q,par);
PhiRS2=bbPhidot2(i,j,uz2pr,s1pr,s2pr,tn,q,par);
Phik=[PhiRS1;PhiRS2];
Phi=Add(Phi,Phik,m,0);
m=m+2;
end

k=k+1;
end

%Euler Parameter Normalization Constraints
i=1;
while i<=nb
[r1,p1]=qPart(q,i);
Phik=(p1'*p1-1)/2;
Phi=Add(Phi,Phik,m,0);
i=i+1;
m=m+1;
end

end




