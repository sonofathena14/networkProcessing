function [arcs, nodes, path]=fix_degree_2(arcs_old, nodes_old, path_old, rootID)

% Find degree 2 nodes to fix
for i=1:length(nodes_old)
    if nodes_old(i,5)==2
        if nodes_old(i,1)==rootID
            continue
        end
        nodeID=nodes_old(i,1);
        vessels=[];
        for j=1:length(arcs_old)
            if arcs_old{1,j}(1,1)==nodeID || arcs_old{1,j}(1,2)==nodeID
                vessels=[vessels; j];
            end
        end
        if length(vessels)~=2
            % Was labeled degree 2 but has more than 2 vessels attached
            error(['Node ', num2str(nodeID), ' was labeled degree 2 but has ' , num2str(length(vessels)), ' vessels attached.'])
        else
            % Truly a degree 2 node
            v1ID=vessels(1);  v2ID=vessels(2);
            [new_vessel]=connectAtNode(arcs_old, v1ID, v2ID, nodeID);
            arcs_old{1, v1ID}=new_vessel;
            arcs_old{1, v2ID}=[];
            node_row=find(nodes_old(:,1)==nodeID);
            nodes_old(node_row, :)=zeros([1, 5]);
            
            for k=1:length(path_old)
                if path_old{1,k}(1,end)==nodeID
                    path_old{1,k}=[];
                elseif ismember(nodeID, path_old{1,k})
                    [~, loc] = ismember(nodeID, path_old{1,k});
                    path_old{1,k}(:,loc)=[];
                end
            end
            path_old=path_old(~cellfun('isempty',path_old));
            
            disp(['Vessels ', num2str(v1ID), ' & ', num2str(v2ID), ' have been merged at node ', num2str(nodeID), '.'])
        end
    end
    arcs_old=arcs_old(~cellfun('isempty',arcs_old));
end

nodes_old( ~any(nodes_old,2), : ) = [];

nodes=nodes_old;
arcs=arcs_old;
path=path_old;

end