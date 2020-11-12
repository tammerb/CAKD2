function A=AT(p)
%eulerToA  Calculates A from euler parameters in p according to Eq 2.5.22
%   [3*3 Orientation Matrix A] = AT(4*1 euler parameters p)
e0=p(1);                    % split e and e0 in p
e=[p(2);p(3);p(4)];
A=(e0^2-e'*e)*eye(3)+2*e*e'+2*e0*atil(e);
end