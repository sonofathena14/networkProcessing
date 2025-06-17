%This function removes blobs in data that result from error in slicing. It
%also serves as a backup to remove duplicated vessels, if missed in
%correctionEngine
%
%   Inputs:
% 1. arcs_old, output from formatData script (or any correction algorithm)
% 2. nodes_old, output from formatData script (or any correction algorithm)
% 3. path_old, output from correctionEngine
%
% Outputs:
% 1. path, {1 x # nodes} cell array, cell path{1, i} is a list of nodes in
%    the shortest path to node i
% 2. arcs, updated arcs cell array
% 3. nodes, updated nodes matrix

% NOTE: This function calls the "dijkstra" function as well as
% "findNearestPoints"

function [path,arcs,nodes]=blobRemover_mice(nodes_old,arcs_old,path_old, rootID,sf)
BR = 0;
segments=zeros(length(arcs_old),3); %extacts VES ID and the input/output nodes
for i=1:length(arcs_old)
    segments(i,1)=i;
    segments(i,2)=arcs_old{1,i}(1,1);
    segments(i,3)=arcs_old{1,i}(1,2);
end


%% Remove blobs
ans1=input('Is the blob at a junction?','s');% different procedures with and without junction
if ans1~='Y'
    node1=input('Input node of blob to remove: '); %Node closest to root
    node2=input('Output node of blob to remove: '); %node on other side of blob
    segments=zeros(length(arcs_old),3); %Redefine vessels
    for i=1:length(arcs_old)
        segments(i,1)=i;
        segments(i,2)=arcs_old{1,i}(1,1);
        segments(i,3)=arcs_old{1,i}(1,2);
    end
    [~, spath] = dijkstra(nodes_old,segments,node1,node2); %use dijkstra to identify the shortest path through blob
    nearest_nodes = findNearestPointsMice(nodes_old,sf); %identify all nodes in blob (within x voxels from center node)
    nodes_to_remove = nearest_nodes(~ismember(nearest_nodes,spath)); %if the node is off path, store to remove it
    vessels_to_remove = [];
    % remove vessels connected to one of the nodes that need to be removed
    for i=1:length(nodes_to_remove) %cycle each node to be removed
        node1 = nodes_to_remove(i);
        for j=1:length(nearest_nodes) %cycle every node
            node2 = nearest_nodes(j);
            if node1 == node2
                continue
            end
            if node1 == rootID || node2 == rootID
                continue
            end
            for x=1:length(segments) %search vessels for 
                if node1==segments(x,2) && node2 == segments(x,3)
                    vessels_to_remove = [vessels_to_remove; x];
                    break
                end
                if node1==segments(x,3) && node2 == segments(x,2)
                    vessels_to_remove = [vessels_to_remove; x];
                    break
                end
            end
        end
    end
    vessels_to_remove = unique(vessels_to_remove,'rows');
    [num_to_remove,~]=size(vessels_to_remove);
    for i=1:num_to_remove %removes vessels and corrects degree of nodes
        vessel_ind=vessels_to_remove(i,1);
        node1=arcs_old{1,vessel_ind}(1,1) ;
        node2=arcs_old{1,vessel_ind}(1,2) ;
        node_ind1=find(nodes_old(:,1)==node1);
        node_ind2=find(nodes_old(:,1)==node2);
        nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
        nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
        arcs_old{1,vessel_ind}=[];
        disp(['Edge from node ', num2str(node1), ' to node ', num2str(node2), ' has been removed.'])
    end
    arcs_old=arcs_old(~cellfun('isempty',arcs_old));
    for i=1:length(nodes_to_remove) %removes necessary nodes
        node_row = find(nodes_old(:,1)==nodes_to_remove(i));
        if nodes_old(node_row,5)~=0
            continue
        end
        nodes_old(node_row, :)=zeros([1, 5]);
        disp(['Node ', num2str(nodes_to_remove(i)), ' has been removed.'])
    end
end
if ans1 =='Y'
    node1=input('Input node of blob to remove: '); %Node closest to root
    nodes = input('[node1 node2] of outputs: ');
    nodes2 = nodes(1); nodes3 = nodes(2);
    [~, spath] = dijkstra(nodes_old,segments,node1,nodes2); %use dijkstra to identify the shortest path through blob in one direction
    [~, path2] = dijkstra(nodes_old,segments,node1,nodes3); %use dijkstra to identify the shortest path through blob in other direction
    nearest_nodes = findNearestPoints(nodes_old,sf); %identify all nodes in blob (within x voxels from center node)
    nodes_to_remove = nearest_nodes(~ismember(nearest_nodes,spath)&~ismember(nearest_nodes,path2)); %if the node is off path, store to remove it
    vessels_to_remove = [];
    % remove vessels connected to one of the nodes that need to be removed
    for i=1:length(nodes_to_remove) %cycle each node to be removed
        node1 = nodes_to_remove(i);
        for j=1:length(nearest_nodes) %cycle every node
            node2 = nearest_nodes(j);
            if node1 == node2
                continue
            end
            if node1 == rootID || node2 == rootID
                continue
            end
            for x=1:length(segments) %search vessels for 
                if node1==segments(x,2) && node2 == segments(x,3)
                    vessels_to_remove = [vessels_to_remove; x];
                    break
                end
                if node1==segments(x,3) && node2 == segments(x,2)
                    vessels_to_remove = [vessels_to_remove; x];
                    break
                end
            end
        end
    end
    vessels_to_remove = unique(vessels_to_remove,'rows');
    [num_to_remove,~]=size(vessels_to_remove);
    for i=1:num_to_remove %removes stored vessels and corrects node degree
        vessel_ind=vessels_to_remove(i,1);
        node1=arcs_old{1,vessel_ind}(1,1) ;
        node2=arcs_old{1,vessel_ind}(1,2) ;
        node_ind1=find(nodes_old(:,1)==node1);
        node_ind2=find(nodes_old(:,1)==node2);
        nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
        nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
        arcs_old{1,vessel_ind}=[];
        disp(['Edge from node ', num2str(node1), ' to node ', num2str(node2), ' has been removed.'])
    end
    arcs_old=arcs_old(~cellfun('isempty',arcs_old));
    for i=1:length(nodes_to_remove)  %remove all store nodes
        node_row = find(nodes_old(:,1)==nodes_to_remove(i));
        if nodes_old(node_row,5)~=0
            continue
        end
        nodes_old(node_row, :)=zeros([1, 5]);
        disp(['Node ', num2str(nodes_to_remove(i)), ' has been removed.'])
    end
end
[arcs_old, nodes_old, path_old]=fix_degree_2(arcs_old, nodes_old, path_old,rootID);

arcs = arcs_old;
nodes = nodes_old;
path = path_old;
end