function K = Keval(apr,b)
K=[apr'*b,apr'*atil(b);atil(apr)*b,apr*b'+b*apr'-apr'*b*eye(3)];


end

