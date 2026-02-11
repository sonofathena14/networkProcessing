% Script to orient network vessel files from Networks or Pruned folders
% Handles nodesC2 (Networks) and nodesC3 (Pruned) automatically
% Aligns pc1 direction to z-axis without flipping

function [] = orient_all_networks(file_paths, root_node_ids, ploton)
% INPUTS:
%   file_paths - cell array of file paths, can be in two formats:
%                Format 1: 'Pruned/Pruned_Network_m1p4_053107' (relative path without .mat)
%                Format 2: 'Networks/Network_Vessels_m1p4_060407' (relative path without .mat)
%   root_node_ids - (optional) array of root node IDs corresponding to each file in file_paths
%                   Used to ensure root node is in positive z direction
%                   Can be empty [] if not needed
%   ploton - (optional) 1 to visualize each alignment, 0 to skip (default: 0)
%
% EXAMPLES:
%   file_paths = {
%       'Pruned/Pruned_Network_m1p4_053107',
%       'Networks/Network_Vessels_m1p4_060407'
%   };
%   root_ids = [12, 14];  % Root node IDs for each file
%   orient_all_networks(file_paths, root_ids, 0);
%
%   % Or without root node checking:
%   orient_all_networks(file_paths, [], 0);

if nargin < 3
    ploton = 0;
end

if nargin < 2
    root_node_ids = [];
end

% Validate input
if ~iscell(file_paths)
    error('file_paths must be a cell array of strings');
end

% Check root_node_ids size if provided
if ~isempty(root_node_ids) && length(root_node_ids) ~= length(file_paths)
    error('root_node_ids array must have same length as file_paths');
end

num_files = numel(file_paths);
fprintf('Processing %d network files...\n\n', num_files);

% Process each file
for i = 1:num_files
    current_path = file_paths{i};
    fprintf('[%d/%d] Processing: %s\n', i, num_files, current_path);
    
    try
        % Add .mat extension if not present
        if ~endsWith(current_path, '.mat')
            input_file = [current_path, '.mat'];
        else
            input_file = current_path;
        end
        
        if ~isfile(input_file)
            warning('File not found: %s. Skipping...', input_file);
            continue;
        end
        
        % Determine folder type and variable name
        if contains(current_path, 'Pruned')
            var_name = 'nodesC3';
            output_prefix = 'Oriented';
        else  % Networks folder
            var_name = 'nodesC2';
            output_prefix = 'Oriented';
        end
        
        % Load ALL variables from the file
        all_data = load(input_file);
        
        % Extract the nodes data
        nodes_data = all_data.(var_name);
        
        % Get root node ID if provided
        root_node_id = [];
        if ~isempty(root_node_ids)
            root_node_id = root_node_ids(i);
        end
        
        % Extract xyz coordinates (columns 2:4, assuming format: ID|x|y|z|degree)
        xyz_original = nodes_data(:, 2:4);
        
        % Apply alignment to z-axis
        [xyz_aligned, flip_needed] = align_points_to_z(xyz_original, nodes_data(:,1), root_node_id, current_path, ploton);
        
        % Update with aligned coordinates
        nodes_data(:, 2:4) = xyz_aligned;
        
        % Prepare output filename
        % Extract the base name and create output path in Oriented folder
        [~, name, ~] = fileparts(input_file);
        
        % Create Oriented folder if it doesn't exist
        if ~exist('Oriented', 'dir')
            mkdir('Oriented');
        end
        
        output_file = fullfile('Oriented', sprintf('%s_%s.mat', output_prefix, extract_base_name(name)));
        
        % Prepare data for saving: add nodesC4, remove old nodesC2/nodesC3
        save_data = rmfield(all_data, var_name);  % Remove old variable
        save_data.nodesC4 = nodes_data;  % Add new variable
        
        % Save all variables
        save(output_file, '-struct', 'save_data');
        
        if flip_needed
            fprintf('  ✓ Saved to: %s (variable: nodesC4) [root-aligned]\n\n', output_file);
        else
            fprintf('  ✓ Saved to: %s (variable: nodesC4)\n\n', output_file);
        end
        
    catch ME
        warning('Error processing %s: %s\n', current_path, ME.message);
        fprintf('  Skipping to next file...\n\n');
    end
end

fprintf('All networks processed!\n');

end


function base_name = extract_base_name(filename)
% Extract base name from different filename formats
% Examples:
%   'Pruned_Network_m1p4_053107' -> 'm1p4_053107'
%   'Network_Vessels_m1p4_060407' -> 'm1p4_060407'

    % Remove common prefixes
    base_name = strrep(filename, 'Pruned_Network_', '');
    base_name = strrep(base_name, 'Network_Vessels_', '');
    
    % If still has underscores at the start, keep cleaning
    while startsWith(base_name, '_')
        base_name = base_name(2:end);
    end
end


function [aligned_points, flip_needed] = align_points_to_z(points, node_ids, root_node_id, network_name, ploton)
% Aligns point cloud's principal direction to z-axis using PCA
% Ensures root node (if provided) is in positive z direction

flip_needed = false;

% Center the point cloud
points_centered = points - mean(points, 1);

% Run PCA
[coeff, ~, ~] = pca(points_centered);

% First principal component (direction of greatest variance)
pc1 = coeff(:, 1);  % A unit vector

% Desired direction (positive z-axis)
z_axis = [0; 0; 1];

% Ensure pc1 points in the positive z direction (no flipping)
if dot(pc1, z_axis) < 0
    pc1 = -pc1;  % Flip pc1 if it points downward
end

% Compute rotation axis (cross product) and angle
v = cross(pc1, z_axis);
s = norm(v);
c = dot(pc1, z_axis);

% Compute rotation matrix
if s == 0 && c == 1
    % Already aligned
    R = eye(3);
elseif s == 0 && c == -1
    % Opposite direction (shouldn't happen after flip check, but handle it)
    % Rotate 180° around x-axis
    R = [1 0 0; 0 -1 0; 0 0 -1];
else
    % Rodrigues rotation matrix
    vx = [   0   -v(3)  v(2);
           v(3)   0   -v(1);
          -v(2) v(1)    0];
    R = eye(3) + vx + vx^2 * ((1 - c) / (s^2));
end

% Apply rotation
aligned_points = (R * points_centered')';

% Check root node position if provided
if ~isempty(root_node_id)
    % Find the root node
    root_idx = find(node_ids == root_node_id, 1);
    
    if isempty(root_idx)
        warning('Root node ID %d not found in network. Skipping root check.', root_node_id);
    else
        % Get root node z-coordinate
        root_z = aligned_points(root_idx, 3);
        
        % If root is in negative z, flip the entire cloud around xy-plane
        if root_z < 0
            aligned_points(:, 3) = -aligned_points(:, 3);
            flip_needed = true;
            fprintf('  Root node was in negative z (%.2f), flipped cloud\n', root_z);
        else
            fprintf('  Root node in positive z (%.2f), no flip needed\n', root_z);
        end
    end
end

% Verify alignment direction
[coeff_check, ~, ~] = pca(aligned_points);
pc1_check = coeff_check(:, 1);
fprintf('  PC1 z-component after alignment: %.4f (target: ~1.0)\n', abs(pc1_check(3)));

% Plot if requested
if ploton
    figure('Name', sprintf('Alignment: %s', network_name), 'Position', [100 100 1200 500]);
    
    subplot(1, 2, 1);
    plot3(points_centered(:,1), points_centered(:,2), points_centered(:,3), 'b.', 'MarkerSize', 3);
    hold on;
    % Draw PC1 vector
    quiver3(0, 0, 0, pc1(1)*50, pc1(2)*50, pc1(3)*50, 'r', 'LineWidth', 2);
    
    % Highlight root node if provided
    if ~isempty(root_node_id)
        root_idx = find(node_ids == root_node_id, 1);
        if ~isempty(root_idx)
            plot3(points_centered(root_idx,1), points_centered(root_idx,2), ...
                  points_centered(root_idx,3), 'go', 'MarkerSize', 12, 'LineWidth', 2);
        end
    end
    
    title('Original (Centered)');
    xlabel('X'); ylabel('Y'); zlabel('Z'); 
    axis equal; grid on; view(3);
    if ~isempty(root_node_id)
        legend('Points', 'PC1 direction', 'Root node');
    else
        legend('Points', 'PC1 direction');
    end
    
    subplot(1, 2, 2);
    plot3(aligned_points(:,1), aligned_points(:,2), aligned_points(:,3), 'r.', 'MarkerSize', 3);
    hold on;
    % Draw z-axis
    quiver3(0, 0, 0, 0, 0, 50, 'g', 'LineWidth', 2);
    
    % Highlight root node if provided
    if ~isempty(root_node_id)
        root_idx = find(node_ids == root_node_id, 1);
        if ~isempty(root_idx)
            plot3(aligned_points(root_idx,1), aligned_points(root_idx,2), ...
                  aligned_points(root_idx,3), 'go', 'MarkerSize', 12, 'LineWidth', 2);
        end
    end
    
    title('Aligned to Z-Axis');
    xlabel('X'); ylabel('Y'); zlabel('Z'); 
    axis equal; grid on; view(3);
    if ~isempty(root_node_id)
        legend('Points', 'Z-axis', 'Root node');
    else
        legend('Points', 'Z-axis');
    end
    
    sgtitle(sprintf('Network: %s', network_name), 'Interpreter', 'none');
end

end