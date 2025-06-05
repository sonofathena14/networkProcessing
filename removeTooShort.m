function [path,arcs,nodes] = removeTooShort(nodes_old,arcs_old,path_old,rootID,sf)
%removeTooShort cycles through all of the edges and removes any edge that
%is shorter than 5 voxels
%
% Inputs:
% 1. arcs_old, output from formatData script (or any correction algorithm)
% 2. nodes_old, output from formatData script (or any correction algorithm)
% 3. path_old, output from correctionEngine
%
% Outputs:
% 1. path, {1 x # nodes} cell array, cell path{1, i} is a list of nodes in
%    the shortest path to node i
% 2. arcs, updated arcs cell array
% 3. nodes, updated nodes matrix

segments=zeros(length(arcs_old),3); %creates a list of ves id with from and to nodes
for i=1:length(arcs_old)
    segments(i,1)=i;
    segments(i,2)=arcs_old{1,i}(1,1);
    segments(i,3)=arcs_old{1,i}(1,2);
end
vessels_to_remove = [];
nodes_to_remove = [];
for i=1:length(segments)%cycles through all vessels
    node1 = segments(i,2);
    node2 = segments(i,3);
    if node1 == rootID || node2 == rootID
        continue
    end
    node1row = find(nodes_old(:,1) == node1);
    node2row = find(nodes_old(:,1) == node2);
    if nodes_old(node1row,5)==1 || nodes_old(node2row,5)==1 %ensures that one vessel is terminal
        x = nodes_old(node1row,2) - nodes_old(node2row,2);
        y = nodes_old(node1row,3) - nodes_old(node2row,3);
        z = nodes_old(node1row,4) - nodes_old(node2row,4);
        if sqrt(abs(x)^2 + abs(y)^2 + abs(z)^2) < 5*sf %if less than 5 voxels, remove vessel
            vessels_to_remove  = [vessels_to_remove; i];
            if nodes_old(node1row,5)==1 %stores the terminal node to remove it
                nodes_to_remove = [nodes_to_remove; node1];
            end
            if nodes_old(node2row,5)==1
                nodes_to_remove = [nodes_to_remove; node2];
            end
        end
    end
end
[num_to_remove,~]=size(vessels_to_remove);
for i=1:num_to_remove %removes the vessels, same as in correction algorithm
    vessel_ind=vessels_to_remove(i,1);
    node1=arcs_old{1,vessel_ind}(1,1) ;
    node2=arcs_old{1,vessel_ind}(1,2) ;
    node_ind1=find(nodes_old(:,1)==node1);
    node_ind2=find(nodes_old(:,1)==node2);
    nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
    nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
    arcs_old{1,vessel_ind}=[];
    disp(['Edge from node ', num2str(node1), ' to node ', num2str(node2), ' has been removed for being too short.'])
end
arcs_old=arcs_old(~cellfun('isempty',arcs_old)); %resets the arcs to remove the blank vessel spaces
for i=1:length(nodes_to_remove) %removes the terminal vessels
    node_row = find(nodes_old(:,1)==nodes_to_remove(i));
    nodes_old(node_row, :)=zeros([1, 5]);
    disp(['Node ', num2str(nodes_to_remove(i)), ' has been removed.'])
end
nodes_old( ~any(nodes_old,2), : ) = []; %resets the nodes to remove blank rows
[arcs_old, nodes_old, path_old]=removeDegree0(arcs_old, nodes_old, path_old,rootID);
[arcs, nodes, path]=fix_degree_2(arcs_old, nodes_old, path_old,rootID); %removes any newly created degree 2 vessels
end