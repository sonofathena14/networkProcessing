% This function plots arcs and nodes data from the formatData program. The
% data is plotted in a 3D point cloud with node IDs labeled. Examining this
% plot alongside the image in 3D slicer helps identify the root node. Also,
% nodes of degree 4 re labeled in red because those are where possible
% loops/trifurcations occur, so they are of interest.

% Inputs:
% 1. arcs, output from formatData script
% 2. nodes, output from formatData script

% Output:
% Figure--3D plot of point cloud

function []=plotSlicerData_prunedOverOrig_2(arcsC3, vessels_to_prune, Data)

map_IDs=Data.mapIDs;

% Create vectors (x, y, z) of point coordinates to plot
x=[];
y=[];
z=[];
% Create vectors (xP, yP, zP) of pruned point coordinates to plot
xP=[];
yP=[];
zP=[];

NumVessels=length(arcsC3);
for i=1:NumVessels
    [numPoints,~]=size(arcsC3{1,i});
    newX=arcsC3{1,i}(2:numPoints,1);
    newY=arcsC3{1,i}(2:numPoints,2);
    newZ=arcsC3{1,i}(2:numPoints,3);
    x=[x;newX];
    y=[y;newY];
    z=[z;newZ];

    % Find the correct adjusted vessel ID to compare to vessels_to_prune list
    [R, ~] = find(map_IDs(:,2)==i);
    adj_ID = map_IDs(R,1);

    if ~ismember(adj_ID, vessels_to_prune)
        xP=[xP;newX];
        yP=[yP;newY];
        zP=[zP;newZ];
    end
end

z=-1*z;
zP=-1*zP;

orig_coords = [x y z];
pruned_coords = [xP yP zP];
orig_coords2plot = setdiff(orig_coords, pruned_coords, 'rows');

% Plot data, labeling nodes

hold on
title('Pruned Over Original')
scatter3(orig_coords2plot(:,1),orig_coords2plot(:,2),orig_coords2plot(:,3),250,'b','.');
scatter3(pruned_coords(:,1),pruned_coords(:,2),pruned_coords(:,3),250,[0,0.8,1],'.');
view(0,0)
axis equal
hold off

end