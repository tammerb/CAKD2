function K = KT(apr,b)
%Evaluate K(apr,b) of Eq.(2.6.37), given apr and b
K=2*[apr'*b,apr'*atil(b);atil(apr)*b,apr*b'+b*apr'-apr'*b*eye(3)];


end

