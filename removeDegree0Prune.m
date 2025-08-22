function [arcs, nodes]=removeDegree0Prune(arcs_old, nodes_old)

% Removes degree 0 nodes
for i=1:length(nodes_old)
    if nodes_old(i,5)==0
        disp(['Node ', num2str(nodes_old(i,1)), ' with degree 0 has been removed.'])
        nodes_old(i, :)=zeros([1, 5]);
    end
end
nodes_old( ~any(nodes_old,2), : ) = [];
arcs = arcs_old;
nodes = nodes_old;
end