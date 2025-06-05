function [bifurcation_angles]=extract_angles_lee(P, D1, D2, A)

if size(P, 1) ~=3
    P=P';
end
if size(D1, 1) ~=3
    D1=D1';
end
if size(D2, 1) ~=3
    D2=D2';
end
if size(A, 1) ~=3
    A=A';
end

% Vessel trajectories from first bifurcation of mouse p4H1
P=P/norm(P);
D1=D1/norm(D1);
D2=D2/norm(D2);
A=A/norm(A);

% Standard reference frame
E1=[1,0,0]';
E2=[0,1,0]';
E3=[0,0,1]';

% Current frame for this junction
e3 = -P/norm(P);
e1 = A-dot(A, e3)*e3;
e1 = e1/norm(e1);
e2 = cross(e3, e1);

% Compute angles
% (i) Angles psi_1 and psi_2
psi_1=acos(dot(e3, D1)/(norm(e3)*norm(D1)));
psi_2=acos(dot(e3, D2)/(norm(e3)*norm(D2)));

% (ii) Angles theta_1 and theta_2
if dot(D1, e2)==0 && dot(D1, e1)==0
    theta_1 = 0;
else
    theta_1=atan(dot(D1, e2)/dot(D1, e1));
end

if dot(D2, e2)==0 && dot(D2, e1)==0
    theta_2 = 0;
else
    theta_2=atan(dot(D2, e2)/dot(D2, e1));
end

R=zeros([3,3]);
R(1,1)=dot(e1,E1);    R(1,2)=dot(e1,E2);    R(1,3)=dot(e1,E3);
R(2,1)=dot(e2,E1);    R(2,2)=dot(e2,E2);    R(2,3)=dot(e2,E3);
R(3,1)=dot(e3,E1);    R(3,2)=dot(e3,E2);    R(3,3)=dot(e3,E3);

D1_std = R * D1;
D2_std = R * D2;

% Shift theta_1 and 2
if D1_std(1) < 0 
    theta_1 = theta_1 + pi;
elseif D1_std(2) < 0
    theta_1 = theta_1 + 2*pi;
end

if D2_std(1) < 0 
    theta_2 = theta_2 + pi;
elseif D2_std(2) < 0
    theta_2 = theta_2 + 2*pi;
end


% Output bifurcation_angles
bifurcation_angles = [round(psi_1, 6), round(psi_2, 6), round(theta_1, 6), round(theta_2, 6)];

% Check for NaN values:
if isnan(psi_1)
    disp(['psi_1 between vectors P=', num2str(P'), ' and D1=', num2str(D1'), ' is NaN.'])
end
if isnan(psi_2)
    disp(['psi_2 between vectors P=', num2str(P'), ' and D2=', num2str(D2'), ' is NaN.'])
end
if isnan(theta_1)
    disp(['theta_1 for vector D1=', num2str(D1'), ' from vector e1=', num2str(e1'), ' about e3=', num2str(e3'), ' is NaN.'])
end
if isnan(theta_2)
    disp(['theta_2 for vector D2=', num2str(D2'), ' from vector e1=', num2str(e1'), ' about e3=', num2str(e3'), ' is NaN.'])
end


end


