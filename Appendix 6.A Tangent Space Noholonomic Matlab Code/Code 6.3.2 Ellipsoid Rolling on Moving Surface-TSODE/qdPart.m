function [rd,pd,aprd,x1d,y1d]=qdPart(qd)
%Fill in all qd components above and below; example is rolling ellipsoid
rd=[qd(1);qd(2);qd(3)];
pd=[qd(4);qd(5);qd(6);qd(7)];
aprd=[qd(8);qd(9);qd(10)];
x1d=qd(11);
y1d=qd(12);

%add other components of  qd
end