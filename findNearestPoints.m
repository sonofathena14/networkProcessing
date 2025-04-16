function result = findNearestPoints(nodes)
%This program identifies the nodes within 5 units of the given node
%   This will pass back the nodes that need to be removed in blobRemover
%   This REQUIRES the Statistics Toolbox in Matlab
    maxDistance = 10; %find all points within this distance, change as needed but better to keep small
% Extract X, Y, Z columns
    targetID = input('Node in center of blob: '); %node that will identify what should be removed
    X = nodes(:, 2);
    Y = nodes(:, 3);
    Z = nodes(:, 4);
    % Create a KD-tree from the coordinates
    points = [X, Y, Z];
    kdtree = KDTreeSearcher(points);
    % Find the row corresponding to the target ID
    targetRow = nodes(nodes(:, 1) == targetID, :);
    if isempty(targetRow)
        error('Target ID not found in the data.');
    end
    % Extract target X, Y, Z coordinates
    targetPoint = targetRow(1, 2:4);
    % Query the KD-tree for all points within maxDistance
    idx = rangesearch(kdtree, targetPoint, maxDistance);
    % Retrieve the points based on indices
    nearbyPoints = nodes(idx{1}, :);  % idx{1} gives the indices of points within range
    % Return the result (all IDs, X, Y, Z within the max distance)
    result = nearbyPoints(:,1);
end