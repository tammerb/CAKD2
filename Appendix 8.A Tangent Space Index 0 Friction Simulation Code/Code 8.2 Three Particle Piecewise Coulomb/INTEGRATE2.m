function[q,qd,qdd,Lam2,v2,v2d,v2dd,u2,q0,U2,V2,B2,jRepar,npar,Rnorm,...
    CondJ,Sw12,jiter,Cr]=...
    INTEGRATE2(n,Q,Qd,Qdd,Vv2,Vv2d,Vv2dd,LLam2,Uu2,q0,U2,V2,B2,...
    par,jRepar,npar,Cr,mode,Sw12,q1b,q2b,q3b)
   
[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

   if Cr>1

% Parameterization

[v2nm,v2dnm,v2ddnm,q0,U2,V2,B2,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar,mode);

u2nm=zeros(nu2,1);

npar=n-1;
jReparrpt(n)=jRepar;
   end
   
Cr=0;
 
% Integration

if n-1==npar
v2dd=v2ddnm;
if Sw12>1
    err=1;
    Lam2=zeros(nu2,1);
    q=Q(:,n-1);
    qd=Qd(:,n-1);
    qddnm=Qdd(:,n-1);
    
    if mode==2
    qdd=[0;qddnm(2);qddnm(3)];
    end
    if mode==3
    qdd=[qddnm(1);0;qddnm(3)];
    end
    if mode==4
    qdd=[qddnm(1);qddnm(2);0];
    end
    
    M=MEval(q,par);
    Phiq=PhiqEval(q,par,mode);
    while err>0.0001
    QA=QAEval(q,qd,Lam2,par,mode);
    Lam2=B2*Phiq*(QA-M*qdd);
    err=norm(M*qdd+Phiq'*Lam2-QA);
    end
    Sw12=0;
end
end

if n-1>npar
    v2nm=Vv2(:,n-1);
    v2dnm=Vv2d(:,n-1);
    v2ddnm=Vv2dd(:,n-1);
    v2dd=v2ddnm;
    Lam2=LLam2(:,n-1);
    u2nm=Uu2(:,n-1);
end

v2d=v2dnm+(h/2)*(v2ddnm+v2dd);
v2=v2nm+h*v2dnm+((h^2)/4)*(v2ddnm+v2dd);

Lam2=LLam2(:,n-1);

%Jacobian Evaluation
[u2,Iteru] = usolv(u2nm,v2,q0,V2,U2,B2,par,mode,q1b,q2b,q3b);
q=q0+V2*v2-U2*u2;
[B2,Biter]=CorrectB(q,B2,U2,par,mode);
Phiq=PhiqEval(q,par,mode);
D2=(eye(nq)-U2*B2*Phiq)*V2;
qd=D2*v2d;
M=MEval(q,par);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam2,par,mode);
[Rvdd,Rvd,Rv,Phiq]=JacobFull(v2,v2d,v2dd,Lam2,u2,par,q0,V2,U2,B2,...
    mode,q1b,q2b,q3b,Sw12);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv,Phiq'-QAsLam];
CondJ=cond(J);

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol
    
%Residual Calculation
[R,B2,h]=Resid(v2dd,Lam2,v2nm,v2dnm,v2ddnm,u2,B2,q0,V2,U2,par,mode,...
    q1b,q2b,q3b);

if i==1
    R1n=norm(R);
end

% Newton Correction

x=-J\R;

delv2dd=x(1);
delLam2=[x(2);x(3)];
v2dd=v2dd+delv2dd;
Lam2=Lam2+delLam2;

%Evaluate v and vd
v2d=v2dnm+(h/2)*(v2ddnm+v2dd);
v2=v2nm+h*v2dnm+((h^2)/4)*(v2ddnm+v2dd);

err=norm(R);
i=i+1;

    if i>6
err=intol/2;
    end

end

jiter=i-1;
Rnorm=norm(R);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evaluate q

u2=Uu2(:,n-1);
[u2,Iteru]=usolv(u2,v2,q0,V2,U2,B2,par,mode,q1b,q2b,q3b);
Iterurpt(n)=Iteru;
q=q0+V2*v2-U2*u2;
qnormrpt(n)=norm(q);

%Update B
[B2,Biter]=CorrectB(q,B2,U2,par,mode);
Biterrpt(n)=Biter;

% Evaluate qd
Phiq = PhiqEval(q,par,mode);
D2=(eye(nq)-U2*B2*Phiq)*V2;
qd=D2*v2d;

% Evaluate qdd

[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par,mode);
qdd=D2*v2dd-U2*B2*Gam;

end


