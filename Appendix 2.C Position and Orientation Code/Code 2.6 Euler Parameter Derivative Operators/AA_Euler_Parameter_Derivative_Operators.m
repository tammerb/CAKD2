
%MATLAB Code Implementing Euler Parameter Derivative Operators of 
%Sections 2.5 and 2.6 is Contained in Functions that can be Copied
%into Application Code

%Example: f(p)=B'(p,a)*A(p)*b; fsp=K(a,(A(p)*b)+B'(p,a)*B(p,b)
a=[1;1;0];
b=[0;1;1];
p=[1;0;0;0];

f=BT(p,a)'*AT(p)*b;
fsp=KT(a,AT(p)*b)+BT(p,a)'*BT(p,b);