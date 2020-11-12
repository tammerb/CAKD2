function K = Keval(apr,b)
K=2*[apr'*b,apr'*atileval(b);atileval(apr)*b,apr*b'+b*apr'-apr'*b*eye(3)];


end

