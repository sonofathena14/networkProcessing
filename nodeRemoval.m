function [arcs,nodes,path] = nodeRemoval(arcs_old,nodes_old,path_old)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% REMOVE NODES

node_to_remove = input('nodeID of node to remove: ');
node_row = find(nodes_old(:,1)==node_to_remove);
nodes_old(node_row, :)=zeros([1, 5]);
disp(['Node ', num2str(node_to_remove), ' has been removed.'])
nodes_old( ~any(nodes_old,2), : ) = [];
[arcs, nodes, path]=fix_degree_2(arcs_old, nodes_old, path_old);
end