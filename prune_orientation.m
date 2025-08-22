function [orientation]=prune_orientation(arcs, nodes)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numVessels=length(arcs);
pathLength=length(path);
orientation=zeros(numVessels, 1);

for i=1:numVessels
    firstNode=arcs{1,i}(1,1);
    firstNRow = find(nodes(:,1)==firstNode);
    sidx = nodes(firstNRow,2);
    sidy = nodes(firstNRow,3);
    sidz = nodes(firstNRow,4);
    if arcs{1,i}(2,1)==sidx && arcs{1,i}(2,2)==sidy && arcs{1,i}(2,3)==sidz
        orientation(i,1)=1;
    else
        orientation(i,1)=-1;
    end
end


end