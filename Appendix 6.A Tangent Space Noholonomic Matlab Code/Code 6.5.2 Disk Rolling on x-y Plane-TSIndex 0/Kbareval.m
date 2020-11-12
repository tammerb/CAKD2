function Kbar=Kbareval(a,b)

Kbar=2*[a'*b,a'*atil(b);atil(a)*b,a*b'+b*a'-a'*b*eye(3)];


end

