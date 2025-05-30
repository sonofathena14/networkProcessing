function [path,arcs,nodes] = generateTrifurcation(nodes_old,arcs_old,path_old,x)

nodes = input('Input nodes to be merged start with node to be removed [node1 node2]: '); 
nodee = input('End node for vessel to be moved to trifurcation: ');

% find xyz coordinate and radii at new trifurcation node
for i=1:length(arcs_old)
   if arcs_old{1,i}(1, 1:2) == nodes 
         xyzr = arcs_old{1,i}(end, 1:4);
   elseif arcs_old{1,i}(1, 1:2)==flip(nodes)
         xyzr = arcs_old{1,i}(2, 1:4);        
   end
end
idxT = find(nodes_old(:,1)==nodes(2));
idxM = find(nodes_old(:,1)==nodes(1));

nodes_old(idxT,5) = 4;
nodes_old(idxM,5) = 2;
for i=1:length(path_old)
    if path_old{1,i}(end) == nodee
        path_old{1,i}(end-1) = nodes(2);
    end
    if path_old{1,i}(end) == nodes(1)
        path_old{1,i} = [];
    end
end
path_old=path_old(~cellfun('isempty',path_old));

% move vessel to trifurcation node
for i=1:length(arcs_old)
   if arcs_old{1,i}(1, 1:2)==[nodes(1) nodee] 
       arcs_old{1,i}(1,1:2)
       arcs_old{1,i}(1,1) = nodes(2);
       arcs_old{1,i}(2, 1:4) = xyzr;
   elseif arcs_old{1,i}(1, 1:2)==flip([nodes(1) nodee])
       arcs_old{1,i}(1,1:2)
       arcs_old{1,i}(1,2) = nodes(2);
       arcs_old{1,i}(end, 1:4) = xyzr;
   end
end

%disp('Remove node before trifurcation')  
%[arcsN,nodesN,pathN]=nodeRemoval(arcs_old,nodes_old,path_old);

[arcs, nodes, path] = fix_degree_2(arcs_old, nodes_old, path_old,x);

end