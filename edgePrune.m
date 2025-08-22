function [arcs,nodes] = edgePrune(arcs_old,nodes_old,node1,node2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Initialize variables
vessels_to_remove = [node1 node2];
%disp(class(arcs_old))
for i=1:length(arcs_old)
    if arcs_old{1,i}(1, 1:2)==vessels_to_remove | arcs_old{1,i}(1, 1:2)==flip(vessels_to_remove)
        arcs_old{1,i}=[];
        node1_row=find(nodes_old(:,1)==node1);
        node2_row=find(nodes_old(:,1)==node2);
        nodes_old(node1_row, 5)=nodes_old(node1_row, 5)-1;
        nodes_old(node2_row, 5)=nodes_old(node2_row, 5)-1;
        disp(['Edge from node ', num2str(node1), ' to node ', num2str(node2), ' has been removed.'])
        break
    end
end
arcs_old=arcs_old(~cellfun('isempty',arcs_old));
% 
% % remove any empty cells.
arcs=arcs_old(~cellfun('isempty',arcs_old));
nodes = nodes_old;
%auto removes the node if degree 0
[arcs, nodes]=removeDegree0Prune(arcs, nodes);
end