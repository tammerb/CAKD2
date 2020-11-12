function[q,qd,qdd,Lam1,v1,v1d,v1dd,u1,q0,U1,V1,B1,jRepar,npar,...
    Rnorm,CondJ,Sw21,jiter,Cr]=...
INTEGRATE1(n,Q,Qd,Qdd,Vv1,Vv1d,Vv1dd,LLam1,Uu1,q0,U1,V1,B1,par,jRepar,...
npar,Cr,mode,q1b,q2b,q3b,Sw21)
   
[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

   if Cr>1

% Parameterization

[v1nm,v1dnm,v1ddnm,q0,U1,V1,B1,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar,mode);

u1nm=zeros(nu1,1);

npar=n-1;
jReparrpt(n)=jRepar;
Cr=0;
   end
     

%Jacobian Calculation

Lam1nm=LLam1(:,n-1);
if Sw21>1
    qnm=Q(:,n-1);
    qdnm=Qd(:,n-1);
    qddnm=Qdd(:,n-1);
    M=MEval(qnm,par);
    QA=QAEval(qnm,qdnm,Lam1nm,par,mode);
    Phiq=PhiqEval(qnm,par,mode);
    Lam1nm=B1*Phiq*(QA-M*qddnm);
    Sw21=0;
end

qdnm=Qd(:,n-1);

% Integration

if n-1==npar
v1dd=v1ddnm;
Lam1=Lam1nm;
end

if n-1>npar 
v1nm=Vv1(:,n-1);
v1dnm=Vv1d(:,n-1);
v1ddnm=Vv1dd(:,n-1);
v1dd=v1ddnm;
Lam1nm=LLam1(:,n-1);
Lam1=Lam1nm;
u1nm=Uu1(:,n-1);
end

[Rvdd,Rvd,Rv,Phiq]=JacobFull(v1nm,v1dnm,v1ddnm,Lam1nm,u1nm,par,...
    q0,V1,U1,B1,mode,q1b,q2b,q3b,Sw21);

v1d=v1dnm+(h/2)*(v1ddnm+v1dd);
v1=v1nm+h*v1dnm+((h^2)/4)*(v1ddnm+v1dd);

%Jacobian Evaluation
[u1,Iteru] = usolv(u1nm,v1,q0,V1,U1,B1,par,mode,q1b,q2b,q3b);
q=q0+V1*v1-U1*u1;
Phiq=PhiqEval(q,par,mode);
D1=(eye(nq)-U1*B1*Phiq)*V1;
qd=D1*v1d;
M=MEval(q,par);
[B1,Biter]=CorrectB(q,B1,U1,par,mode);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam1,par,mode);
[Rvdd,Rvd,Rv,Phiq]=JacobFull(v1,v1d,v1dd,Lam1,u1,par,...
    q0,V1,U1,B1,mode,q2b,q3b,Sw21);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv,Phiq'-QAsLam];
CondJ=cond(J);

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol;

    
%Residual Calculation
[R,B1,h]=Resid(v1dd,Lam1,v1nm,v1dnm,v1ddnm,u1,B1,q0,V1,U1,par,...
    mode,q1b,q2b,q3b);

if i==1
    R1n=norm(R);
end

% Newton Correction

x=-J\R;

delv1dd=[x(1);x(2)];
delLam1=x(3);
v1dd=v1dd+delv1dd;
Lam1=Lam1+delLam1;

%Evaluate v and vd
v1d=v1dnm+(h/2)*(v1ddnm+v1dd);
v1=v1nm+h*v1dnm+((h^2)/4)*(v1ddnm+v1dd);

err=norm(R);
i=i+1;

    if i>6
err=intol/2;
    end

end

jiter=i-1;
Rnorm=norm(R);

%Evaluate q

[u1,Iteru]=usolv(u1nm,v1,q0,V1,U1,B1,par,mode,q1b,q2b,q3b);
Iterurpt(n)=Iteru;
q=q0+V1*v1-U1*u1;
qnormrpt(n)=norm(q);

%Update B
[B1,Biter]=CorrectB(q,B1,U1,par,mode);
Biterrpt(n)=Biter;

% Evaluate qd
Phiq = PhiqEval(q,par,mode);
D1=(eye(nq)-U1*B1*Phiq)*V1;
qd=D1*v1d;

% Evaluate qdd

[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par,mode);
qdd=D1*v1dd-U1*B1*Gam;


end

