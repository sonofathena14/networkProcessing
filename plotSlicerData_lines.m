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

function []=plotSlicerData_lines(arcs, nodes, color,fignb)

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
figure(fignb)
hold on
scatter3(x,y,z,500,color,'.');
t=text(xN,yN,zN,label,'FontSize',10);
for i=1:size(nodes, 1)
%     if ismember(nodes(i,1),B(:,1))
%         t(i).Color = 'red';
%         t(i).FontSize = 24;
%         t(i).FontWeight = 'bold';
      if nodes(i, 5)==1
        t(i).Color = 'black';
        t(i).FontSize = 18;
      elseif nodes(i, 5)==3
        t(i).Color = 'red';
        t(i).FontSize = 18;
      else
        t(i).Color = 'magenta';
        t(i).FontSize = 24;
      end
%     if nodes(i,1)==302
%         t(i).Color = 'magenta';
%         t(i).FontSize = 18;
%         t(i).FontWeight = 'bold';        
%     elseif nodes(i,1)==303
%         t(i).Color = 'magenta';
%         t(i).FontSize = 18;
%         t(i).FontWeight = 'bold';        
%     end     
end
for i=1:length(arcs)
    A=arcs{1,i};
    [rows, ~]=size(A);
    for j=3:rows
        v1=A(j-1,1:3);
        v2=A(j,1:3);
        v=[v2;v1];
        plot3(v(:,1),v(:,2),v(:,3),'Color', [0 0.4470 0.7410]);
    end
end
%  LABEL COORDINATES
% for i=1:size(x, 1)
%     coordLabel=['(' num2str(x(i)) ', ' num2str(y(i)) ', ' num2str(z(i)) ')'];
%     t2=text(x(i),y(i),z(i),coordLabel,'FontSize',10);
%     if mod(i,2)==0
%         t2.Color = 'black';
%     else
%         t2.Color = 'red';
%     end
%     t2.FontSize = 12;
% end
hold off

% figure(332)
% hold on
% s=scatter3(x,y,z,500,'.','MarkerEdgeColor',[.5 .5 .5]);
% % for i=1:length(arcs)
% %     A=arcs{1,i};
% %     [rows, ~]=size(A);
% %     for j=3:rows
% %         v1=A(j-1,1:3);
% %         v2=A(j,1:3);
% %         v=[v2;v1];
% %         plot3(v(:,1),v(:,2),v(:,3),'r');
% %     end
% % end
% hold off
end