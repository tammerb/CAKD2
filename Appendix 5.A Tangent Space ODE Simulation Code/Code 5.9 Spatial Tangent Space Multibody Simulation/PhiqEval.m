function Phiq=PhiqEval(tn,q,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

Phiq=zeros(nc,ngc);

I3=eye(3);
k=1;        %Joint No.
m=0;        %Constraint equation Counter-1
while k<=nh

%Distance Constraint
if SJDT(1,k)==1
[i,j,s1pr,s2pr,d]=DistPart(k,SJDT);

[Phiq1,Phiq2]=bbPhiqdist(i,j,s1pr,s2pr,d,tn,q,par);
Phiq=Add(Phiq,Phiq1,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiq2,m,7*(j-1));
end

m=m+1;
end

%Spherical Constraint
if SJDT(1,k)==2
[i,j,s1pr,s2pr]=SphPart(k,SJDT);

[Phiq1,Phiq2]=bbPhiqsph(i,j,s1pr,s2pr,tn,q,par);
Phiq=Add(Phiq,Phiq1,m,7*(i-1));      

if j>=1
Phiq=Add(Phiq,Phiq2,m,7*(j-1));
end

m=m+3;
end

%Cylindrical Constraint
if SJDT(1,k)==3
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=CylPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);

Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41];
Phiq2k=[Phiqcyl12;Phiqcyl22;Phiqcyl32;Phiqcyl42];

Phiq=Add(Phiq,Phiq1k,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiq2k,m,7*(j-1));
end

m=m+4;
end

%Revolute Constraint
if SJDT(1,k)==4
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=RevPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
[Phiqrev51,Phiqrev52]=bbPhiqdot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);

Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41;Phiqrev51];
Phiq2k=[Phiqcyl12;Phiqcyl22;Phiqcyl32;Phiqcyl42;Phiqrev52];

Phiq=Add(Phiq,Phiq1k,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiq2k,m,7*(j-1));
end

m=m+5;
end

%Translational Constraint
if SJDT(1,k)==5
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=TranPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;

[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
[Phiqtran51,Phiqtran52]=bbPhiqdot1(i,j,vy1pr,vx2pr,tn,q,par);

Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41;Phiqtran51];
Phiq2k=[Phiqcyl12;Phiqcyl22;Phiqcyl32;Phiqcyl42;Phiqtran52];

Phiq=Add(Phiq,Phiq1k,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiq2k,m,7*(j-1));
end

m=m+5;
end

%Universal Constraint
if SJDT(1,k)==6
[i,j,s1pr,s2pr,vz1pr,vz2pr]=UnivPart(k,SJDT);    
[Phiqks1,Phiqks2]=bbPhiqsph(i,j,s1pr,s2pr,tn,q,par);
[Phiqkd11,Phiqkd12]=bbPhiqdot1(i,j,vz1pr,vz2pr,tn,q,par);
Phiqk1=[Phiqks1;Phiqkd11];
Phiqk2=[Phiqks2;Phiqkd12];
Phiq=Add(Phiq,Phiqk1,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiqk2,m,7*(j-1));
end
m=m+4;
end

%Strut Constraint
if SJDT(1,k)==7
[i,j,s1pr,s2pr,vx1pr,vy1pr]=StrutPart(k,SJDT);
[Phiqstr11,Phiqstr12]=bbPhiqdot2(i,j,vx1pr,s1pr,s2pr,tn,q,par);
[Phiqstr21,Phiqstr22]=bbPhiqdot2(i,j,vy1pr,s1pr,s2pr,tn,q,par);
Phiqk1=[Phiqstr11;Phiqstr21];
Phiqk2=[Phiqstr12;Phiqstr22];
Phiq=Add(Phiq,Phiqk1,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiqk2,m,7*(j-1));
end
m=m+2;
end

%Revolute-Spherical Constraint
if SJDT(1,k)==8
[i,j,s1pr,s2pr,d,vz2pr]=RevSphPart(k,SJDT);
[PhiqRS11,PhiqRS12]=bbPhiqdist(i,j,s1pr,s2pr,d,tn,q,par);
[PhiqRS21,PhiqRS22]=bbPhiqdot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);
Phiqk1=[PhiqRS11;PhiqRS21];
Phiqk2=[PhiqRS12;PhiqRS22];
Phiq=Add(Phiq,Phiqk1,m,7*(i-1));

if j>=1
Phiq=Add(Phiq,Phiqk2,m,7*(j-1));
end
m=m+2;
end

k=k+1;

end

%Euler Parameter Normalization Constraints
i=1;
while i<=nb
[r1,p1]=qPart(q,i);
Phiq1=[zeros(1,3),p1'];
Phiq=Add(Phiq,Phiq1,m,7*(i-1));
m=m+1;
i=i+1;
end

end


