function [arcs] = OrientArcs(arcs,nodes)

disp('in OrientArcs')
arcsT = arcs;
for i = 1:length(arcs)
    N = length(arcs{1,i}(:,1));
    if abs(sum(arcs{1,i}(2,1:3)-arcs{1,i}(3,1:3))) > 3
       i
       arcsT{1,i}(2,:) = arcs{1,i}(N,:);
       arcsT{1,i}(N,:) = arcs{1,i}(2,:);
       arcsT{1,i}(1,1) = arcs{1,i}(1,2);
       arcsT{1,i}(1,2) = arcs{1,i}(1,1);
       arcs = arcsT;
    end
end
