function [V,SSD, xyzl]=linReg(xyz, hubNode_coordinates, plot_on)

xyz=xyz';
hubNode_coordinates=hubNode_coordinates';

% Engine  
% xyz0 = mean(xyz,2);
xyz0=hubNode_coordinates;
if xyz(:, 1) ~= xyz0
    xyz = fliplr(xyz);
end

if size(xyz, 2) >= 3
    A = xyz-xyz0;
    [U,S,~] = svd(A);
    V = U(:,1);
    SSD=S(2,2)^2+S(3,3)^2;
elseif size(xyz, 2) == 2
    V = xyz(:, 1) - xyz(:, 2);
    SSD = 0;
else
    error('Cannot fit a line to just one point.')
end

t = [0 1];
t1 = min(t);
t2 = max(t);
xyzl = xyz0 + [t1,t2].*V; % size 3x2

% Make sure vector is properly aligned with points
test_vector = xyz(:,end) - xyz(:,1);
DP = dot(V, test_vector);
check = 1 - DP;
if check >=1
    % Vector is misaligned
    V = -V;
end

if plot_on == 1
    % % Check
    x = xyz(1,:);
    y = xyz(2,:);
    z = xyz(3,:);
    xl = xyzl(1,:);
    yl = xyzl(2,:);
    zl = xyzl(3,:);
    %close all
    hold on
    plot3(x,y,z,'ob');
    plot3(xl,yl,zl,'r','LineWidth', 2);
%     scatter3(xl(1),yl(1),zl(1),'*m','LineWidth', 5);
%     scatter3(xl(2),yl(2),zl(2),'*c','LineWidth', 5);
end

end