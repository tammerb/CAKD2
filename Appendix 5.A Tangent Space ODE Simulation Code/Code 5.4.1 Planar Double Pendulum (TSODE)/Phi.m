function Phi = Phi(t,q,par)

Phi=[q(1)-cos(q(3));q(2)-sin(q(3));q(1)+cos(q(3))-q(4)+cos(q(6));...
    q(2)+sin(q(3))-q(5)+sin(q(6))];
end

