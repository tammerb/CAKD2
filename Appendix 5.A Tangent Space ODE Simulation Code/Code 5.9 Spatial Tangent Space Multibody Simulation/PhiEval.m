function Phi=PhiEval(tn,q,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

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
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=CylPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

Phicyl1=bbPhidot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
Phicyl2=bbPhidot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
Phicyl3=bbPhidot1(i,j,vz1pr,vx2pr,tn,q,par);
Phicyl4=bbPhidot1(i,j,vz1pr,vy2pr,tn,q,par);

Phik=[Phicyl1;Phicyl2;Phicyl3;Phicyl4];
Phi=Add(Phi,Phik,m,0);
m=m+4;
end

%Revolute Constraint
if SJDT(1,k)==4
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=RevPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

Phicyl1=bbPhidot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
Phicyl2=bbPhidot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
Phicyl3=bbPhidot1(i,j,vz1pr,vx2pr,tn,q,par);
Phicyl4=bbPhidot1(i,j,vz1pr,vy2pr,tn,q,par);
Phirev5=bbPhidot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);

Phik=[Phicyl1;Phicyl2;Phicyl3;Phicyl4;Phirev5];
Phi=Add(Phi,Phik,m,0);
m=m+5;
end

%Translational Constraint
if SJDT(1,k)==5
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=TranPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

Phicyl1=bbPhidot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
Phicyl2=bbPhidot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
Phicyl3=bbPhidot1(i,j,vz1pr,vx2pr,tn,q,par);
Phicyl4=bbPhidot1(i,j,vz1pr,vy2pr,tn,q,par);
Phitran5=bbPhidot1(i,j,vy1pr,vx2pr,tn,q,par);

Phik=[Phicyl1;Phicyl2;Phicyl3;Phicyl4;Phitran5];
Phi=Add(Phi,Phik,m,0);
m=m+5;
end

%Universal Constraint
if SJDT(1,k)==6
[i,j,s1pr,s2pr,vz1pr,vz2pr]=UnivPart(k,SJDT);    
Phiks=bbPhisph(i,j,s1pr,s2pr,tn,q,par);
Phikd1=bbPhidot1(i,j,vz1pr,vz2pr,tn,q,par);
Phik=[Phiks;Phikd1];
Phi=Add(Phi,Phik,m,0);
m=m+4;
end

%Strut Constraint
if SJDT(1,k)==7
[i,j,s1pr,s2pr,vx2pr,vy2pr]=StrutPart(k,SJDT);
Phistr1=bbPhidot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
Phistr2=bbPhidot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
Phik=[Phistr1;Phistr2];
Phi=Add(Phi,Phik,m,0);
m=m+2;
end

%Revolute-Spherical Constraint
if SJDT(1,k)==8
[i,j,s1pr,s2pr,d,vz2pr]=RevSphPart(k,SJDT);
PhiRS1=bbPhidist(i,j,s1pr,s2pr,d,tn,q,par);
PhiRS2=bbPhidot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);
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




