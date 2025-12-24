% This algorithm performs point cloud  alignment using ICP initialization
% and CPD refinement. It preserves original sizes (NO SCALING), NO PCA
% pre-alignment.

function [aligned_cloud,tform, alignment_error] = align_vascular_icp_cpd_no_pca_mod(second_cloud, aligned_cloud_arg, ploton)

%% Read in the point cloud to be aligned

folder = 'Point_clouds_inferior';

M2_file = fullfile(folder, strcat(second_cloud, '_points_cloud.csv'));

M2 = readmatrix(M2_file);

%% Center point cloud

M2_centered = M2 - mean(M2, 1);

%% Stage 1: ICP initialization (no PCA)

fixed_pc  = pointCloud(aligned_cloud_arg);
moving_pc = pointCloud(M2_centered);

try
    tform_icp = pcregistericp(moving_pc, fixed_pc, ...
        'Metric','pointToPlane', ...
        'MaxIterations',30, ...
        'Tolerance',[0.001, 0.001]);

    M2_icp_aligned = pctransform(moving_pc, tform_icp).Location;
    fprintf('ICP initialization completed\n');
catch ME
    fprintf('ICP failed: %s\nUsing raw centered cloud\n', ME.message);
    M2_icp_aligned = M2_centered;
end

%% Stage 2: CPD refinement

moving_pc = pointCloud(M2_icp_aligned);

cpd_options = struct(...
    'Transform', 'rigid', ...     % no scaling
    'MaxIterations', 150, ...
    'Tolerance', 1e-5, ...
    'OutlierRatio', 0.2, ...
    'Verbose', true ...
);

try
    tform = pcregistercpd(moving_pc, fixed_pc, cpd_options);
    aligned_pc = pctransform(moving_pc, tform);
    aligned_cloud = aligned_pc.Location;
    fprintf('CPD refinement completed successfully\n');
catch ME
    fprintf('CPD failed: %s\nUsing ICP alignment only\n', ME.message);
    aligned_cloud = M2_icp_aligned;
    tform = [];
end

%% Stage 3: Alignment quality metric (no subsampling)

alignment_error = calculate_alignment_error(aligned_cloud_arg, aligned_cloud);
fprintf('Final alignment error: %.4f\n', alignment_error);

% Save results
writematrix(aligned_cloud, 'alignedM2_icp_cpd.csv');

%% Visualization

if ploton == 1
    figure('Position', [100, 100, 1200, 400]);

    % Original centered clouds
    subplot(1,3,1);
    plot3(aligned_cloud_arg(:,1), aligned_cloud_arg(:,2), aligned_cloud_arg(:,3), 'b.', 'MarkerSize', 6);
    title('M1 Centered');
    axis equal; grid on;

    subplot(1,3,2);
    plot3(M2_centered(:,1), M2_centered(:,2), M2_centered(:,3), 'r.', 'MarkerSize', 6);
    title('M2 Centered');
    axis equal; grid on;

    % Final alignment
    subplot(1,3,3);
    plot3(aligned_cloud_arg(:,1), aligned_cloud_arg(:,2), aligned_cloud_arg(:,3), 'b.', 'MarkerSize', 6);
    hold on;
    plot3(aligned_cloud(:,1), aligned_cloud(:,2), aligned_cloud(:,3), 'r.', 'MarkerSize', 6);
    title(sprintf('Final Alignment (Error: %.4f)', alignment_error));
    legend('M1','M2 (Final)');
    axis equal; grid on;

    sgtitle('Vascular Network Alignment (ICP + CPD)', 'FontSize', 16);
end

end 

function error = calculate_alignment_error(fixed_points, moving_points)
% Calculate alignment quality using bidirectional nearest neighbor distance
    D1 = pdist2(fixed_points, moving_points, 'euclidean','Smallest',1);
    D2 = pdist2(moving_points, fixed_points, 'euclidean','Smallest',1);
    error = (mean(D1) + mean(D2)) / 2;
end
