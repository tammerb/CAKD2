function K=KEval(apr,b)
K=2*[apr'*b,apr'*atil(b);atil(apr)*b,apr*b'+b*apr'-apr'*b*eye(3)];


end

