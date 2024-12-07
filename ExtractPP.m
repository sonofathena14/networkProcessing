% This program extracts the principal pathway from a network.
% The principal pathway is defined to contain the largest connected subtree
% of vessels with radii above a certain percentage of the MPA radius.

function [arcs_pp, nodes_pp, path_pp, newNetwork_pp, Connection_pp,...
        arcs_offpp, nodes_offpp, newNetwork_offpp, Connection_offpp]=ExtractPP(arcs, nodes, newNetwork, MPA_ID)

% Initialize variables for principal pathway and off principal pathway data.
arcs_pp=cell(0);
nodes_pp=[];
arcs_offpp=cell(0);
nodes_offpp=[];
newNetwork_offpp=[];

% Obtain information about MPA
rootRad=newNetwork(MPA_ID,4);
rootStartNode=newNetwork(MPA_ID, 2);

% You can reset the ratio threshold here. Currently at 40%.
ratio_threshold = 0.4;

% Cycle through newNetwork of original tree
for i=1:length(newNetwork)
    rtr=newNetwork(i,4)/rootRad;
    if rtr >= ratio_threshold
        % Include vessel in PP
        arcs_pp=[arcs_pp arcs{1,i}];
        topNode=newNetwork(i,2);
        botNode=newNetwork(i,3);
        N=[topNode; botNode];
        
        % C will include the nodes that are new from this vessel
        if ~isempty(nodes_pp)
            C=setdiff(N, nodes_pp(:,1));
        else 
            C=N;
        end
        
        % Include row in nodes_pp for each node in C
        for j=1:length(C)
            id=C(j);
            ind=find(nodes(:,1)==id);
            nodes_pp=[nodes_pp; nodes(ind,:)];
        end
    else
        % Include vessel in offPP
        arcs_offpp=[arcs_offpp arcs{1,i}];
        newNetwork_offpp = [newNetwork_offpp; newNetwork(i,:)]; 
        topNode=newNetwork(i,2);
        botNode=newNetwork(i,3);
        N=[topNode; botNode];
        
        % C will include the nodes that are new from this vessel
        if ~isempty(nodes_offpp)
            C=setdiff(N, nodes_offpp(:,1));
        else 
            C=N;
        end
        
        % Include row in nodes_offpp for each node in C
        for j=1:length(C)
            id=C(j);
            ind=find(nodes(:,1)==id);
            nodes_offpp=[nodes_offpp; nodes(ind,:)];
        end
    end
end

[path_pp]=useDijkstra(arcs_pp, nodes_pp, rootStartNode);

% At this point, principal pathway may have disconnected components. 
% Need to move these to offpp. We detect these by looking for NaN entries in path.
for i=1:length(path_pp)
    if isnan(path_pp{1,i})
        % Node from row i in nodes is not connected
        path_pp{1,i} = [];
        nodeID = nodes_pp(i, 1);
        
        % Add nodeID to nodes_offpp if not already included
        C = setdiff(nodeID, nodes_pp(:, 1));
        if ~isempty(C)
            nodes_offpp=[nodes_offpp; nodes_pp(:, 1)];
        end
        
        % Zero out row for node in nodes_pp
        nodes_pp(i, :) = zeros(1, 5);
        
        % look for any vessels from pp connected to this node.
        % move them to offpp.
        for j=1:length(arcs_pp)
            if ~isempty(arcs_pp{1,j})
                if arcs_pp{1,j}(1,1)==nodeID || arcs_pp{1,j}(1,2)==nodeID
                    
                    % Move arc to arcs_offpp
                    arcs_offpp = [arcs_offpp arcs_pp{1,j}];
                    nodes_currentVessel = arcs_pp{1,j}(1, 1:2);
                    
                    % Move newNetwork_row for this arc to newNetwork_offpp
                    newNetwork_row = find(newNetwork(:, 2:3) == nodes_currentVessel);
                    if isempty(newNetwork_row)
                        newNetwork_row = find(newNetwork(:, 2:3) == fliplr(nodes_currentVessel));
                        newNetwork_offpp = [newNetwork_offpp; newNetwork(newNetwork_row, :)];
                    end
                    
                    % Remove from arcs_pp
                    arcs_pp{1,j} = [];
                end
            end
        end
    end
end

% Remove all zero rows from nodes_pp 
nodes_pp( ~any(nodes_pp,2), : ) = [];

% Remove all empty cells from path_pp and arcs_pp
path_pp = path_pp(~cellfun('isempty',path_pp));
arcs_pp = arcs_pp(~cellfun('isempty',arcs_pp));

% Update degree column of nodes_pp
for i=1:length(nodes_pp)
    nodeID = nodes_pp(i, 1);
    degree = 0;
    
    % Cycle through arcs_pp, adding 1 to degree for each vessel you find.
    for j=1:length(arcs_pp)
        if arcs_pp{1,j}(1,1)==nodeID || arcs_pp{1,j}(1,2)==nodeID
            degree = degree + 1;
        end
    end
    
    % Update degree value in nodes_pp
    nodes_pp(i, 5) = degree;
end

% Update degree column of nodes_offpp
for i=1:length(nodes_offpp)
    nodeID = nodes_offpp(i, 1);
    degree = 0;
    
    % Cycle through arcs_offpp, adding 1 to degree for each vessel you find.
    for j=1:length(arcs_offpp)
        if arcs_offpp{1,j}(1,1)==nodeID || arcs_offpp{1,j}(1,2)==nodeID
            degree = degree + 1;
        end
    end
    
    % Update degree value in nodes_offpp
    nodes_offpp(i, 5) = degree;
end

% Update vessel IDsn in newNetwork_offpp
newNetwork_offpp(:, 1) = 1:1:length(newNetwork_offpp);

[arcs_pp, nodes_pp, path_pp]=fix_degree_2(arcs_pp, nodes_pp, path_pp);
plotSlicerData(arcs_pp, nodes_pp,'b')

[orientation_pp]=edge_orientation(arcs_pp, path_pp);
[newNetwork_pp]=networkGenerator(arcs_pp, orientation_pp);
[Connection_pp, ~, ~]=connectivity(nodes_pp, newNetwork_pp);
[Connection_offpp, ~, ~]=connectivity(nodes_offpp, newNetwork_offpp);
end