% Code that attempts to align a 3D point cloud's direction of greatest variance to the z-axis via
% PCA and Rodrigues rotation algorithm. Currently being used to preprocess point clouds for convex
% hull creation and comparison. 

function [aligned_cloud] = align_cloud_to_z(cloud_file, ploton)

%{
% Load point cloud
filename = strcat(cloud_file, '_points_cloud.csv');
%}
% Fixed folder name (subfolder of Desktop)
folder1 = 'Point_clouds_inferior';
folder2 = 'inferior';

% Read in the full reference cloud 

filename1 = fullfile(folder1, strcat(cloud_file, '_points_cloud.csv'));

% Read in the principal pathway of the reference cloud 

filename2 = fullfile(folder2, strcat(cloud_file, '_principal_pathways.csv'));

%Debugging for file location

if ~isfile(filename1)
    error('File not found: %s\nCurrent folder: %s', filename1, pwd);
end

if ~isfile(filename2)
    error('File not found: %s\nCurrent folder: %s', filename2, pwd);
end

M_principal = readmatrix(filename2);

M_principal(:, end) = [];

M_main = readmatrix(filename1);

% Center full cloud
M_centered = M_main - mean(M_main, 1);

% Center principal pathway cloud
M_principal_centered = M_principal - mean(M_principal, 1);

% Run PCA
[coeff, ~, ~] = pca(M_principal_centered);

% First principal component (direction of greatest variance)
pc1 = coeff(:, 1);  % A unit vector

% Desired direction (z-axis)
z_axis = [0; 0; 1];

% Compute rotation axis (cross product) and angle
v = cross(pc1, z_axis);
s = norm(v);
c = dot(pc1, z_axis);

% If already aligned, no rotation needed
if s == 0 && c == 1
    R = eye(3);
elseif s == 0 && c == -1
    % Opposite direction – rotate 180° around any perpendicular axis
    % Find orthogonal axis (e.g., x or y not colinear with pc1)
    perp = null(pc1'); 
    rot_axis = perp(:,1); 
    R = axis_angle_to_matrix(rot_axis, pi);
else
    %Rodrigues rotation matrix
    vx = [   0   -v(3)  v(2);
           v(3)   0   -v(1);
          -v(2) v(1)    0];
    R = eye(3) + vx + vx^2 * ((1 - c) / (s^2));
end

% Apply rotation
aligned_cloud = (R * M_centered')';

% Save
writematrix(aligned_cloud, 'M1_aligned_to_z.csv');

% Plot
if ploton
    figure;
    subplot(1, 2, 1);
    plot3(M_centered(:,1), M_centered(:,2), M_centered(:,3), 'b.');
    title('Original Centered Cloud');
    xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; grid on;

    subplot(1, 2, 2);
    plot3(aligned_cloud(:,1), aligned_cloud(:,2), aligned_cloud(:,3), 'r.');
    title('Aligned to Z-Axis');
    xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; grid on;
end

end

function R = axis_angle_to_matrix(rot_axis, theta)
    rot_axis = rot_axis / norm(rot_axis);
    x = rot_axis(1); y = rot_axis(2); z = rot_axis(3);
    c = cos(theta); s = sin(theta); C = 1 - c;

    R = [x*x*C + c,   x*y*C - z*s, x*z*C + y*s;
         y*x*C + z*s, y*y*C + c,   y*z*C - x*s;
         z*x*C - y*s, z*y*C + x*s, z*z*C + c];
end
