function [rd,pd,ad,sd]=qdPart(qd)
%Fill in all qd components above and below; example is rolling ellipsoid
rd=[qd(1);qd(2);qd(3)];
pd=[qd(4);qd(5);qd(6);qd(7)];
ad=[qd(8);qd(9)];
sd=qd(10);


%add other components of  qd
end