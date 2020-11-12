function A=Add(A,B,m,n)

%Add elements of matrix B to elements of A below and to the right of (m,n)
[r,s]=size(B);
i=1;
while i<=r
j=1;
while j<=s
A(m+i,n+j)=A(m+i,n+j)+B(i,j);
j=j+1;
end
i=i+1;
end

end


