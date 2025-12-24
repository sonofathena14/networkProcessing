% This code utilizes two separate algorithms to align multiple point clouds
% to a reference point cloud. The reference point cloud is first aligned
% to the z-axis. Then, the user specified point clouds are aligned to the
% z-aligned reference point cloud. The user can then specify whether to
% plot the superimposed clouds to visually check for alignment.

function [] = multiple_cloud_alignment1(ploton)

% --- Reference cloud ---
ref_name = input('Please enter a reference cloud (ex. m1p1_053007): ','s');

% --- List of clouds to align ---
clouds = input('Enter clouds to align (comma-separated): ','s');
cloud_entries = strtrim(strsplit(clouds, ','));

% Align reference to z-axis
ref_cloud_aligned = align_cloud_to_z(ref_name,0);

% Cell array to store aligned clouds
aligned_clouds = cell(numel(cloud_entries),1);

% Align each cloud
for i = 1:numel(cloud_entries)
    currentName = strtrim(cloud_entries{i});
    [aligned_cloud,~,~] = ...
        align_vascular_icp_cpd_no_pca_mod(currentName, ref_cloud_aligned, 0);
    aligned_clouds{i} = aligned_cloud;
end

% --- Plot all clouds superimposed ---
if ploton
    figure('Name','Superimposed Aligned Clouds');
    hold on; grid on; axis equal
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3)

    % Plot reference cloud
    h = gobjects(numel(cloud_entries)+1,1);
    legendNames = cell(numel(cloud_entries)+1,1);

    h(1) = plot3(ref_cloud_aligned(:,1), ...
                 ref_cloud_aligned(:,2), ...
                 ref_cloud_aligned(:,3), ...
                 '.', 'MarkerSize',6, 'Color',[0 0 0]);  % black
    legendNames{1} = sprintf('Reference (%s)',ref_name);

    % Plot each aligned cloud with different color
    cmap = lines(numel(cloud_entries)); % nice distinct colors
    for i = 1:numel(cloud_entries)
        h(i+1) = plot3(aligned_clouds{i}(:,1), ...
                       aligned_clouds{i}(:,2), ...
                       aligned_clouds{i}(:,3), ...
                       '.', 'MarkerSize',6, 'Color',cmap(i,:));
        legendNames{i+1} = strtrim(cloud_entries{i});
    end

    legend(h, legendNames, 'Interpreter','none','Location','bestoutside');
    title('Reference + Aligned Clouds (ICP + CPD)')
end
end
