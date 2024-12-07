% This function plots arcs and nodes data from the formatData program. The
% data is plotted in a 3D point cloud with node IDs labeled. Examining this
% plot alongside the image in 3D slicer helps identify the root node. Also,
% nodes of degree 4 are labeled in red because those are where possible
% loops/trifurcations occur, so they are of interest.
% 
% Inputs:
% 1. arcs, output from formatData script
% 2. nodes, output from formatData script
% 
% Output:
% Figure--3D plot of point cloud

function []=plotSlicerData(arcs, nodes, color,fignb)

% Create vectors (x, y, z) of point coordinates to plot
x=[];
y=[];
z=[];
NumVessels=length(arcs);
for i=1:NumVessels
    [numPoints,~]=size(arcs{1,i});
    newX=arcs{1,i}(2:numPoints,1);
    newY=arcs{1,i}(2:numPoints,2);
    newZ=arcs{1,i}(2:numPoints,3);
    x=[x;newX];
    y=[y;newY];
    z=[z;newZ];
end

% Create vector of node ID labels
label=num2str(nodes(:,1));
xN=nodes(:,2);
yN=nodes(:,3);
zN=nodes(:,4);

% Plot data, labeling nodes
figure(fignb); hold on
scatter3(x,y,z,10,color,'o');

t=text(xN,yN,zN,label,'FontSize',16);
set(gca,'fontsize',16);

for i=1:length(nodes)
    if nodes(i,5)==4
        t(i).Color = 'red';
        t(i).FontSize = 16;
        t(i).FontWeight = 'bold';
    else
        t(i).Color = 'black';
        t(i).FontSize = 16;
    end
end
view(0,0)
hold off
end

% function [] = plotSlicerData(arcs, nodes, color)
% 
% % Create vectors (x, y, z) of point coordinates to plot
% x = [];
% y = [];
% z = [];
% NumVessels = length(arcs);
% 
% figure(1); hold on
% % Plot each vessel as a line
% for i = 1:NumVessels
%     [numPoints, ~] = size(arcs{1, i});
%     newX = arcs{1, i}(2:numPoints, 1);
%     newY = arcs{1, i}(2:numPoints, 2);
%     newZ = arcs{1, i}(2:numPoints, 3);
%     plot3(newX, newY, newZ, 'Color', color, 'LineWidth', 2);
% end
% 
% % Create vector of node ID labels
% label = num2str(nodes(:, 1));
% xN = nodes(:, 2);
% yN = nodes(:, 3);
% zN = nodes(:, 4);
% 
% % Label nodes
% t = text(xN, yN, zN, label, 'FontSize', 10);
% set(gca, 'fontsize', 16);
% 
% for i = 1:length(nodes)
%     if nodes(i, 5) == 4
%         t(i).Color = 'red';
%         t(i).FontSize = 18;
%         t(i).FontWeight = 'bold';
%     else
%         t(i).Color = 'black';
%         t(i).FontSize = 18;
%     end
% end
% 
% view(0, 0)
% hold off
% end
