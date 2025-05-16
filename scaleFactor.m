function [arcs,nodes] = scaleFactor(arcs_old,nodes_old)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
root_vessel=input('[node1 node2] of root edge: ');
for i=1:length(arcs_old)
    if arcs_old{1,i}(1, 1:2)==root_vessel | arcs_old{1,i}(1, 1:2)==flip(root_vessel)
        radii = [];
        for j=1:5
            radii = [radii; arcs_old{1,i}(j,4)];
        end
        ravg = mean(radii);
        scaler = 430/ravg;
        break
    end
end
for i=1:length(arcs_old)
    arcs_old{1,i}(2:end, 1:4)=arcs_old{1,i}(2:end, 1:4) *scaler;
end
nodes_old(:,2:4) = nodes_old(:,2:4)*scaler;

arcs = arcs_old;
nodes = nodes_old;


end