% This program generates the Vessel nodes and arcs for a structured tree
% defined by alpha, beta, lrr, and angles

function [Vessel_arcs, Vessel_nodes]=st_vessels_rmin(rmin, root_start_xyz, root_end_xyz, ...
                                     root_rad, alpha, beta, lrr1, lrr2, Angles)

tic

psi1 = Angles(1); 
psi2 = Angles(2); 
theta1 = Angles(3); 
theta2 = Angles(4);

% Set up arrays Vessel_arcs and Vessel_nodes
Vessel_nodes = [0 root_start_xyz 1; 1 root_end_xyz 3];
length_root = norm(root_end_xyz-root_start_xyz);
Vessel_arcs = [0 1 root_rad length_root];
branching_nodes_visited = [];
branching_nodes_unvisited = 1;
NodeID = 1;

done = 0;
while ~isempty(branching_nodes_unvisited)
    
    % Iterate unvisited branching nodes. Try to create 2 branches off each.
    for i=1:size(branching_nodes_unvisited, 1)
        parent_node = branching_nodes_unvisited(i);
        
        % Parent radius:
        rP = Vessel_arcs(find(Vessel_arcs(:, 2)==parent_node), 3);
        
        % Daughter radii:
        rD1 = alpha * rP;
        rD2 = beta  * rP;
        
        if rD2 < rmin
            % If rD2 is under rmin, don't generate branches.
            % Change degree to degree 1.
            Vessel_nodes( find(Vessel_nodes(:, 1) == parent_node), 5)=1;
            
        else
            
            lrr_d1 = lrr1*exp(-lrr2*rD1);
            lrr_d2 = lrr1*exp(-lrr2*rD2);
            
            % Lengths of daughters:
            lD1 = rD1 * lrr_d1;
            lD2 = rD2 * lrr_d2;
            
            % Parent vector:
            grandparent_node = Vessel_arcs(find(Vessel_arcs(:, 2)==parent_node), 1);
            P_start_xyz = Vessel_nodes(find(Vessel_nodes(:, 1)==grandparent_node), 2:4);
            P_end_xyz   = Vessel_nodes(find(Vessel_nodes(:, 1)==parent_node), 2:4);
            
            P_start_xyz = P_start_xyz';
            P_end_xyz   = P_end_xyz';
            
            P = P_end_xyz-P_start_xyz;
            e3 = P/norm(P);
            
            % Aunt vector:
            grandparents_kids = Vessel_arcs(find(Vessel_arcs(:, 1)==grandparent_node), 2);
            aunt_node = setdiff(grandparents_kids, parent_node);
            A_end_xyz = Vessel_nodes(find(Vessel_nodes(:, 1)==aunt_node), 2:4);
            A_end_xyz = A_end_xyz';
            A  = A_end_xyz-P_start_xyz;
            A  = A/norm(A);
            
            % If no aunt, use [1,0,0]
            if isempty(A)
                A = [1,0,0]';
            end
            
            e1 = A-dot(A, e3)*e3;
            e1 = e1/norm(e1);
            e2 = cross(e3, e1);          

            pick = rand(1);
%             if pick >= 0.12
            if pick >= 0
                theta2 = -theta2;
            end
            
            % Generate in-plane rotations D1_tilde and D2_tilde
            D1_tilde = rotation_matrix(e2,psi1) * e3 * lD1 + P_end_xyz;
            D2_tilde = rotation_matrix(e2,psi2) * e3 * lD2 + P_end_xyz;
            
            % Generate out-of-plane rotations D1 and D2
            D1 = rotation_matrix(e3,theta1) * (D1_tilde - P_end_xyz) + P_end_xyz;
            D2 = rotation_matrix(e3,theta2) * (D2_tilde - P_end_xyz) + P_end_xyz;
            
            % Add nodes and coords to Vessel_nodes matrix and edges to Vessel_arcs 
            % Add Daughter 1
            NodeID = NodeID + 1;
            Vessel_nodes = [Vessel_nodes; NodeID D1' 3];
            Vessel_arcs = [Vessel_arcs; parent_node NodeID rD1 lD1];
            
            % Add Daughter 2
            NodeID = NodeID + 1;
            Vessel_nodes = [Vessel_nodes; NodeID D2' 3];
            Vessel_arcs = [Vessel_arcs; parent_node NodeID rD2 lD2];
            
            % Add parent to branching_nodes_visited
            branching_nodes_visited = [branching_nodes_visited; parent_node];
            
        end
    end
    
    % All nodes with degree=3 are branching nodes.
    branching_nodes = Vessel_nodes( find(Vessel_nodes(:, 5) == 3), :);
    
    % All branching nodes whose branches aren't completed are unvisited.
    branching_nodes_unvisited = setdiff(branching_nodes(:,1), branching_nodes_visited);
    
end

% Add labels for Arcs array
Vessel_arcs = mat2cell(Vessel_arcs, ones(1, size(Vessel_arcs,1)), ones(1, 4));
Vessel_arcs = [cell(1, 4); Vessel_arcs];
Vessel_arcs{1,1} = 'StartNode';
Vessel_arcs{1,2} = 'EndNode';
Vessel_arcs{1,3} = 'Radius';
Vessel_arcs{1,4} = 'Length';

% Add labels for Nodes array
Vessel_nodes = mat2cell(Vessel_nodes, ones(1, size(Vessel_arcs,1)), ones(1, 5));
Vessel_nodes = [cell(1, 5); Vessel_nodes];
Vessel_nodes{1,1} = 'NodeID';
Vessel_nodes{1,2} = 'x';
Vessel_nodes{1,3} = 'y';
Vessel_nodes{1,4} = 'z';
Vessel_nodes{1,5} = 'degree';

plot_3D_tree(Vessel_arcs, Vessel_nodes)

toc

end