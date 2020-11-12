function P2 = P2(t,q,x,par)

[r1,phi1,r2,phi2]=qPart(q);

P2=[0,0,x(3)*cos(phi1),0,0,0;...
    0,0,x(3)*sin(phi1),0,0,0;...
    0,0,-x(3)*cos(phi1),0,0,-x(6)*cos(phi2);...
    0,0,-x(3)*sin(phi1),0,0,-x(6)*sin(phi2)];
end

