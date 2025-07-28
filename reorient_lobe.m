% Code that attempts to align two 3D point clouds along a similar axis via
% PCA and SVD. Currently being used to preprocess point clouds for convex
% hull creation and comparison. 

function [aligned_cloud] = reorient_lobe(first_cloud, second_cloud, ploton)

% Read in the two clouds 

%M1 = strcat(first_cloud, '_convex_hull.csv');

%M2 = strcat(second_cloud, '_convex_hull.csv');

M1 = strcat(first_cloud, '_points_cloud.csv');

M2 = strcat(second_cloud, '_points_cloud.csv');

M1 = readmatrix(M1);

M2 = readmatrix(M2);

% Center both point clouds
M1_centered = M1 - mean(M1, 1);
M2_centered = M2 - mean(M2, 1);

% Then run PCA on centered clouds to find principle directions
[coeff1, ~, ~] = pca(M1_centered);
[coeff2, ~, ~] = pca(M2_centered);

[U, ~, V] = svd(coeff2 * coeff1');

R = V * U'; % Rotation matrix

%Check for reflections

if det(R) < 0
    V(:,3) = -V(:,3);
    R = V * U';
end

% Rotate the second point cloud via the rotation matrix

M2_aligned = (R * M2_centered')'

aligned_cloud = M2_aligned;
%{
% Find principle directions using PCA

[coeff1,~,~] = pca(M1);

[coeff2,~,~] = pca(M2);
%}

% Get rotational components and calculate rotation matrix (Kabsch
% algorithm)

[U, ~, V] = svd(coeff2 * coeff1');

R = V * U'; % Rotation matrix

%Check for reflections

if det(R) < 0
    V(:,3) = -V(:,3);
    R = V * U';
end

% Rotate the second point cloud via the rotation matrix

M2_aligned = (R * M2_centered')'

aligned_cloud = M2_aligned;

% Save aligned point cloud

writematrix(M2_aligned, 'alignedM2.csv');

% Plot point clouds for visual inspection

%{
if ploton == 1
    figure;
    %plot3(M1(:,1), M1(:,2), M1(:,3), 'b.')
    plot3(M1_centered(:,1), M1_centered(:,2), M1_centered(:,3), 'b.')

    hold on
    plot3(M2_aligned(:,1), M2_aligned(:,2), M2_aligned(:,3), 'r.')
    legend('First Point Cloud', 'Aligned Second Cloud','FontSize', 14);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on; axis equal;
    title('Point Cloud Alignment', 'FontSize', 18);
end

%}

if ploton == 1
    figure;

    % Subplot 1: First centered point cloud
    subplot(1,3,1); % 1 row, 3 columns, 1st plot
    plot3(M1_centered(:,1), M1_centered(:,2), M1_centered(:,3), 'b.')
    title('M1 Centered');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;

    % Subplot 2: Second centered point cloud
    subplot(1,3,2); % 2nd plot
    plot3(M2_aligned(:,1), M2_aligned(:,2), M2_aligned(:,3), 'r.')
    title('M2 Centered');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;

    % Subplot 3: Superimposed original M1 and aligned M2
    subplot(1,3,3); % 3rd plot
    plot3(M1_centered(:,1), M1_centered(:,2), M1_centered(:,3), 'b.')
    hold on
    plot3(M2_aligned(:,1), M2_aligned(:,2), M2_aligned(:,3), 'r.')
    title('M1 vs Aligned M2');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('M1', 'M2 Aligned');
    axis equal; grid on;
end



    


    

