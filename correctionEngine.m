% This program removes vessels from the original network so that no false
% branches, duplicate edges, duplicate points, or cycles exist in the final tree

% Inputs:
% 1. arcs_old, output from formatData script
% 2. nodes_old, output from formatData script
% 3. root_nodeID, identified by examining the plot

% Outputs:
% 1. path, {1 x # nodes} cell array, cell path{1, i} is a list of nodes in
%    the shortest path to node i
% 2. arcs, updated arcs cell array
% 3. nodes, updated nodes matrix
% 4. correction_log, table of errors of each type

% NOTE: This function calls the "dijkstra" function.

function [path, arcs, nodes, correction_log]=correctionEngine(arcs_old,nodes_old,root_nodeID)

FB=0;
DE=0;
DP=0;
SC=0;

%% DOUBLE EDGES
dup_vessels=[];
dup_indices=[];

for i=1:length(arcs_old)
    edge=arcs_old{1,i}(1,1:2);  % For each edge...
    for j=1:length(arcs_old)
        if j~=i & arcs_old{1,j}(1,1:2)==edge % ...if another entry contains the same edge,
            if isempty(dup_vessels)==1
                % First duplicate found
                dup_vessels=[dup_vessels; edge];
                dup_indices=[dup_indices; i j];
            elseif ~isempty(dup_vessels) & ismember(edge,dup_vessels,'rows')==0
                % Subsequent duplicates found
                dup_vessels=[dup_vessels;edge];
                dup_indices=[dup_indices; i j];
                break
            end
        end
    end
end
% Remove the false duplicate, ie: whichever one is longer.
[numDup,~]=size(dup_vessels);
for i=1:numDup
    node1=dup_vessels(i,1);
    vessel_ind1=dup_indices(i,1);
    [size1,~]=size(arcs_old{1,vessel_ind1});
    node2=dup_vessels(i,2);
    vessel_ind2=dup_indices(i,2);
    [size2,~]=size(arcs_old{1,vessel_ind2});
    if size1>=size2
        arcs_old{1,vessel_ind1}=[];
        node_ind1=find(nodes_old(:,1)==node1);
        nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
        node_ind2=find(nodes_old(:,1)==node2);
        nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
        disp(['Double edge from ', num2str(node1), ' to ', num2str(node2), ' has been removed.'])
        DE=DE+1;
    elseif size2>size1
        arcs_old{1,vessel_ind2}=[];
        node_ind1=find(nodes_old(:,1)==node1);
        nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
        node_ind2=find(nodes_old(:,1)==node2);
        nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
        disp(['Double edge from ', num2str(node1), ' to ', num2str(node2), ' has been removed.'])
        DE=DE+1;
    end
end
arcs_old=arcs_old(~cellfun('isempty',arcs_old));

%% DOUBLE EDGESMore actions
segments=zeros(length(arcs_old),3);
for i=1:length(arcs_old)
    segments(i,1)=i;
    segments(i,2)=arcs_old{1,i}(1,1);
    segments(i,3)=arcs_old{1,i}(1,2);
end

duplicate_vessels = [];

% Initialize a set to store already seen routes (both same and reverse)
seen_routes = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Loop through each vessel in the segments table
for i = 1:size(segments, 1)
    vessel_id = segments(i, 1);
    from_node = segments(i, 2);
    to_node = segments(i, 3);

    % Create a key for the same route (from-to)
    route_key = sprintf('%d-%d', from_node, to_node);

    % Create a key for the reverse route (to-from)
    reverse_route_key = sprintf('%d-%d', to_node, from_node);

    % Check if the same route or reverse route has already been seen
    if isKey(seen_routes, route_key) || isKey(seen_routes, reverse_route_key)
        % If it has been seen before, add the current vessel to duplicates list
        duplicate_vessels = [duplicate_vessels; vessel_id];
    else
        % Otherwise, mark this route (both same and reverse) as seen
        seen_routes(route_key) = true;
        seen_routes(reverse_route_key) = true;
    end
end

%Remove vessels
[num_duplicate_vessels,~]=size(duplicate_vessels);
for i=1:num_duplicate_vessels
    vessel_ind=duplicate_vessels(i,1);
    node1=arcs_old{1,vessel_ind}(1,1) ;
    node2=arcs_old{1,vessel_ind}(1,2) ;
    node_ind1=find(nodes_old(:,1)==node1);
    node_ind2=find(nodes_old(:,1)==node2);
    nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
    nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
    arcs_old{1,vessel_ind}=[];
    disp(['Double edge from node ', num2str(node1), ' to node ', num2str(node2), ' has been removed.'])
end
arcs_old=arcs_old(~cellfun('isempty',arcs_old));
nodes_old( ~any(nodes_old,2), : ) = [];
%duplicate edges now removed

%% DUPLICATE POINTS
for i=1:length(arcs_old)
    [n,~]=size(arcs_old{1,i});
    for j=3:n
        currentPoint=arcs_old{1,i}(j,1:3);
        prevPoint=arcs_old{1,i}(j-1,1:3);
        if currentPoint==prevPoint
            arcs_old{1,i}(j,:)=[0 0 0 0];
            disp(['Duplicate edge point in arc ', num2str(i), ' has been removed.'])
            DP=DP+1;
        end
    end
    arcs_old{1,i}( ~any(arcs_old{1,i},2), : ) = []; 
end

%% SHORT CYCLES
% Turn arcs into segments, so they are the appropriate format for inputting
% into the dijkstra algorithm
segments=zeros(length(arcs_old),3);
for i=1:length(arcs_old)
    segments(i,1)=i;
    segments(i,2)=arcs_old{1,i}(1,1);
    segments(i,3)=arcs_old{1,i}(1,2);
end

vessels_to_remove=[];

[numNodes,~]=size(nodes_old);
before=[];
after=[];

% Run Dijkstra to break short cycles. Info on dijkstra algorithm is
% included in the code comments of the function.
[~,path_old]=dijkstra(nodes_old,segments,root_nodeID);

% Check if tree is disconnected. The "path" cell will have an entry of NaN
% if any nodes could not be reached from the rootnode, indicating that our
% tree would be disconnected. Display an error message if this is the case.
nodes_that_cant_be_reached=[];
for i=1:numNodes
    if isnan(cell2mat(path_old(i)))==1 % If path entry is NaN, then that node couldn't be reached
        nodes_that_cant_be_reached=[nodes_that_cant_be_reached; nodes_old(i,1)];
    end
end
if isempty(nodes_that_cant_be_reached)==0
    %nodes_that_cant_be_reached
    msg='The above nodes could not be reached from the rootNode.';
    disp(msg) % Throw a message if any nodes could not be reached from the root node
end

% Now, search through all edges in the arcs structure and remove any edges
% that are not a part of the "path" found by the dijkstra algorithm.
for i=1:length(arcs_old)                  % for each edge...
    edge=arcs_old{1,i}(1,1:2);  % ... edge=[node1 node2].
    if edge(1)~=root_nodeID
        for j=1:numNodes+1
            if j<numNodes+1         % Search through all cells of "path".
                [~,n]=size(path_old{1,j});
                node1_ind=find(path_old{1,j}==edge(1));  % Look for cells containing node1.
                if isempty(node1_ind)==0             % isempty=false, meaning you've found a cell containing node1
                    before=path_old{1,j}(1,node1_ind-1); % before=node listed before node1
                    if before==edge(2)               % If that's node2, you've found the edge in path, so it should be included in the final network.
                        break
                    elseif node1_ind~=n
                        after=path_old{1,j}(1,node1_ind+1); % after=node listed after node1
                        if after==edge(2)               % If that's node2, you've found the edge in path, so it should be included in the final network.
                            break
                        end
                    end
                end
            end
        end
    end    
    if j==numNodes+1  % If you've searched through every cell in "path" the loop hasn't broken, that means the edge wasn't found in the path.
        vessels_to_remove=[vessels_to_remove; i]; % So that edge is part of a cycle and should be removed.
    end
end

% Remove edges excluded by Dijkstra
[num_to_remove,~]=size(vessels_to_remove);
for i=1:num_to_remove
    vessel_ind=vessels_to_remove(i,1);
    node1=arcs_old{1,vessel_ind}(1,1) ;
    node2=arcs_old{1,vessel_ind}(1,2) ;
    node_ind1=find(nodes_old(:,1)==node1);
    node_ind2=find(nodes_old(:,1)==node2);
    nodes_old(node_ind1,5)=nodes_old(node_ind1,5)-1;
    nodes_old(node_ind2,5)=nodes_old(node_ind2,5)-1;
    arcs_old{1,vessel_ind}=[];
    disp(['Small cycle edge from node ', num2str(node1), ' to node ', num2str(node2), ' has been removed.'])
    SC=SC+1;
end
arcs_old=arcs_old(~cellfun('isempty',arcs_old));

%% Fix degree 2 nodes, and collect correction details
[arcs, nodes, path]=fix_degree_2(arcs_old, nodes_old, path_old,root_nodeID);
correction_log = table(FB, DE, DP, SC);



end




