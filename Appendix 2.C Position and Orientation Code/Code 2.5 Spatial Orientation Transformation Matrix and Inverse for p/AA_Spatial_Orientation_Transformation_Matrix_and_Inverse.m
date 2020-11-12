%AA_Spatial_Orientation_Transformation_Matrix_and_Inverse for p

Part=[1];    %Part=1, Define orthogonal tramsformation matrix 
            %Part=2, Compute p for given orthogonal tramsformation matrix A


if Part==1  %Define orthogonal tramsformation matrix

mode=[1];     %mode=1, enter A; mode=2, point definition of A;
             %mode=3, Euler'S Theorem

if mode==1
A=eye(3);
end

if mode==2
%Define orthogonal transformation matrix A using point definition of 
%Section 2.5.5.

rO=[;;];    %Enter vector rO to origin
rP=[;;];    %Enter vector rP to point P
rQ=[;;];    %Enter vector rQ to point Q

f=(1/norm(rP-rO))*(rP-rO);
h=(1/norm(atil(f)*(rQ-rO)))*atil(f)*(rQ-rO);
g=-atil(f)*h;

A=[f,g,h];

end 

if mode==3  %Define orthogonal transformation matrix A using Euler's
            %Theorem of Section 2.5.3.

%Unit vector u about which rotation angle chi brings x-y-z frame into 
%coincidence with x'-y'-z' frame

ubar=[;;];    %Enter vector ubar about which rotation occurs, 
               %not necessarily normalized
               
chi=[];    %Enter angle chi of rotation 

u=(1/norm(ubar))*ubar;
e0=cos(chi/2);
e=sin(chi/2)*u;
A=(e0^2-e'*e)*eye(3)+2*e*e'+e0*atil(e);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Part==2  %Given orthogonal A, compute Euler parameter vector p in 
            %Section 2.5.4 so A=A(p)

[a11,a12,a13,a21,a22,a23,a31,a32,a33] = partA(A);

if norm(A'-A)>0
trA=a11+a22+a33;
e0=0.5*sqrt(trA+1);
e=(1/(2*sqrt(trA+1)))*[a32-a23;a13-a31;a21-a12];
end 

if norm(A'-A)==0
e0=1;
e=zeros(3,1);
end
        
p=[e0;e];
 
end




