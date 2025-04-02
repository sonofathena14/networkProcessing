function R = rotation_matrix(u,theta)
R = zeros(3,3);
ux = u(1); uy = u(2); uz = u(3);
c = cos(theta);
s = sin(theta);

R(1,1) = c+ux.^2 .*(1-c);
R(1,2) = ux.*uy.*(1-c)+uz.*s;
R(1,3) = ux.*uz.*(1-c)-uy.*s;
R(2,1) = ux.*uy.*(1-c)-uz.*s;
R(2,2) = c+uy.^2.*(1-c);
R(2,3) = uy.*uz.*(1-c)+ux.*s;
R(3,1) = ux.*uz.*(1-c)+uy.*s;
R(3,2) = uy.*uz.*(1-c)-ux.*s;
R(3,3) = c + uz.^2.*(1-c);
end