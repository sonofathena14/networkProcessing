function [] = plot_3D_tree(Vessel_arcs, Vessel_nodes)

tic
figure(1)
hold on
grid on
Vessel_nodes = Vessel_nodes([2:end],:);
Vessel_nodes = cell2mat(Vessel_nodes);
maxrad = max(cell2mat(Vessel_arcs([2:end],3)));

maxZ = max(abs(Vessel_nodes(:, 4)));

for i=2:size(Vessel_arcs, 1)
    
    StartNode = Vessel_arcs{i,1};
    EndNode = Vessel_arcs{i,2};
    rad = Vessel_arcs{i,3};
    
    Start_xyz = Vessel_nodes( find(Vessel_nodes(:,1)==StartNode), 2:4);
    End_xyz = Vessel_nodes( find(Vessel_nodes(:,1)==EndNode), 2:4);
    
    width = 7*(rad/maxrad);
    C = [1-Start_xyz(3)/maxZ,0,Start_xyz(3)/maxZ]
    plot3([Start_xyz(1),End_xyz(1)], [Start_xyz(2),End_xyz(2)], [Start_xyz(3),End_xyz(3)], '-', 'LineWidth', width, 'Color', C);


end

view(0,0)

end