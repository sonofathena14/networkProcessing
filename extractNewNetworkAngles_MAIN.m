function [Angles, CP_Groups, CP_Counts, CP_matrix, vessel_vectors] = extractNewNetworkAngles_MAIN(arcs, nodes, newNetwork, conn)

% Create array 'Angles'
Angles = cell(2, 5);
Angles{1,1} = 'hubNode';
Angles{1,2} = 'psi_1';
Angles{1,3} = 'psi_2';
Angles{1,4} = 'theta_1';
Angles{1,5} = 'theta_2';

% Create array 'CP_Groups'
CP_Groups = cell(2, 4);
CP_Groups{1,1} = 'hubNode';
CP_Groups{1,2} = 'Parent_Change_Points';
CP_Groups{1,3} = 'D1_Change_Points';
CP_Groups{1,4} = 'D2_Change_Points';
P_CPs = [];
D1_CPs = [];
D2_CPs = [];

% Create array 'CP_Counts'
CP_Counts = cell(2, 3);
CP_Counts{1,1} = '0_Change_Points';
CP_Counts{1,2} = '1_Change_Points';
CP_Counts{1,3} = '2_Change_Points';

n0cp = 0;   n1cp = 0;   n2cp = 0;

% Create CP_matrix
CP_matrix = cell(4,4);
CP_matrix{2, 1} = 'P';
CP_matrix{3, 1} = 'D1';
CP_matrix{4, 1} = 'D2';
CP_matrix{1, 2} = '0_CP';
CP_matrix{1, 3} = '1_CP';
CP_matrix{1, 4} = '2_CP';

% Create array 'vessel_vectors'
vessel_vectors = cell(2, 5);
vessel_vectors{1,1} = 'hubNode';
vessel_vectors{1,2} = 'P';
vessel_vectors{1,3} = 'D1';
vessel_vectors{1,4} = 'D2';
vessel_vectors{1,5} = 'A';

P_vectors = [];
D1_vectors = [];
D2_vectors = [];
A_vectors = [];

% Create vectors to hold angle values
hubs = [];
psi_1 = [];
psi_2 = [];
theta_1 = [];
theta_2 = [];

% maxD = max # daughters at a single bif in network
maxD = max(conn(:, 2));

for i=1:length(conn)
    daughters=conn(i, 4:3+maxD);
    daughters=daughters(daughters ~= 0);
    numD=length(daughters);
    hubNode=conn(i,1);
    if numD==2 % If we are at a bifurcation:
        
        % Compute angles
        
        hubNode_ind = find(nodes(:, 1) == hubNode);
        hubNode_coordinates = nodes(hubNode_ind, 2:4);
        
        % First, get the three vessels involved.
        % Parent vessel ID
        pID=conn(i,end);
        if pID~=0
            pArc=arcs{1,pID};

            % Daughters 1 & 2 vessel IDs
            d1ID=conn(i,4);
            d2ID=conn(i,5);

            % Make sure D1 has larger radius.
            radD1=newNetwork(d1ID,4);
            radD2=newNetwork(d2ID,4);
            if radD2>radD1
                t=d1ID;
                d1ID=d2ID;
                d2ID=t;
            end

            % OK, now D1 is the daughter of larger radius
            d1Arc=arcs{1,d1ID};
            d2Arc=arcs{1,d2ID};
          
            px = pArc(2:end, 1);
            d1x = d1Arc(2:end, 1);
            d2x = d2Arc(2:end, 1);

            py = pArc(2:end, 2);
            d1y = d1Arc(2:end, 2);
            d2y = d2Arc(2:end, 2);

            pz = pArc(2:end, 3);
            d1z = d1Arc(2:end, 3);
            d2z = d2Arc(2:end, 3);

          
            %Calculating distance of vessel

            xyz_p = [px, py, pz];
            pDistanceCalc = sqrt(sum(diff(xyz_p).^2, 2));  % Euclidean distance between successive points
            pDistance = [0; cumsum(pDistanceCalc)];

            xyz_d1 = [d1x, d1y, d1z];
            d1DistanceCalc = sqrt(sum(diff(xyz_d1).^2, 2));  % Euclidean distance between successive points
            d1Distance = [0; cumsum(d1DistanceCalc)];

            xyz_d2 = [d2x, d2y, d2z];
            d2DistanceCalc = sqrt(sum(diff(xyz_d2).^2, 2));  % Euclidean distance between successive points
            d2Distance = [0; cumsum(d2DistanceCalc)];
            %pause;

    


% 3D Vessel Plots


%{
           
            figure(1);
            scatter3(px, py, pz, 'filled', 'MarkerFaceColor', 'b')
            hold on
            scatter3(d1x, d1y, d1z, 'filled', 'MarkerFaceColor', 'r')
            scatter3(d2x, d2y, d2z, 'filled', 'MarkerFaceColor', 'g')
            legend('pArc', 'd1Arc', 'd2Arc')
            pause;
           
%}


            %figure(1);



% Radius as a function of Euclidean distance plots
    
        %{
            figure(2); clf;
            subplot(3,1,1)
               plot(pDistance,pArc(2:end,4),'o-','linewidth',3,'markersize',7)
               set(gca,'fontsize',16)
               xlabel('distance')
               ylabel('r Parent') 
            subplot(3,1,2)
               plot(d1Distance, d1Arc(2:end,4),'o-','linewidth',3,'markersize',7)
               set(gca,'fontsize',16)
               xlabel('distance')
               ylabel('r Daugther 1')
            subplot(3,1,3)
               plot(d2Distance, d2Arc(2:end,4),'o-','linewidth',3,'markersize',7)
               set(gca,'fontsize',16)
               xlabel('distance')
               ylabel('r Daugther 2')
            pause;
        
        %}

            % Then, get the coordinates from each edge 
            % Parent:
            xyz_P = pArc(2:end, 1:3);
            [V_P, min_AIC, min_BIC, numCPs_P] = chooseVesselVector(xyz_P, hubNode_coordinates);
            
            % Daughter 1:
            xyz_D1 = d1Arc(2:end, 1:3);
            [V_D1, min_AIC, min_BIC, numCPs_D1] = chooseVesselVector(xyz_D1, hubNode_coordinates);

            % Daughter 2:
            xyz_D2 = d2Arc(2:end, 1:3);
            [V_D2, min_AIC, min_BIC, numCPs_D2] = chooseVesselVector(xyz_D2, hubNode_coordinates);    

            % Resize vectors (3 x 1)
            if size(V_P, 1) ~=3
                V_P=V_P';
            end
            if size(V_D1, 1) ~=3
                V_D1=V_D1';
            end
            if size(V_D2, 1) ~=3
                V_D2=V_D2';
            end
            
            % Find Aunt vector
            topNode = newNetwork(pID, 2);
            aunt_conn_row = find(conn(:, 1)==topNode);
            daughters=conn(aunt_conn_row, 4:3+maxD);
            daughters=daughters(daughters ~= 0);
            numD=length(daughters);
            
            % If no Aunt exists:
            if numD == 1
                % Use perpendicular vector to parent in same y-plane

                A = [V_P(3); V_P(2); -V_P(1)];
            else
                daughters = setdiff(daughters, pID);
                aunt_vector = daughters;
                if size(daughters) >= 1
                    lengths = [];
                    for j=1:length(daughters)
                        L = newNetwork(daughters(j), 5);
                        lengths = [lengths; L];
                    end
                    [~, ind] = max(lengths);
                    aunt_vector = daughters(ind);
                end
                auntTop = newNetwork(aunt_vector, 2);
                auntBot = newNetwork(aunt_vector, 3);
                topCoords = nodes( find(nodes(:,1)==auntTop), 2:4);
                botCoords = nodes( find(nodes(:,1)==auntBot), 2:4);
                A = botCoords - topCoords;
                A = A/norm(A);
            end
            
            % Resize A
            if size(A, 1) ~=3
                A=A';
            end
            
            
            % Check if your A works, i.e. does e1 exist how it is
            % currently defined?
            test_e3 = -V_P/norm(V_P);
            A = A/norm(A);
            test_e1 = A - (dot(A, test_e3)*test_e3);
            if round(test_e1, 8)==0


                % Try using perpendicular vector to parent in same x-plane
                A = [V_P(1); V_P(3); -V_P(2)];

                % Check if this works, i.e. does e1 exist how it is
                % currently defined?
                A = A/norm(A);
                test_e1 = A - (dot(A, test_e3)*test_e3);
                if round(test_e1, 8)==0


                    % Try using perpendicular vector to parent in same z-plane
                    A = [V_P(2); -V_P(1); V_P(3)];

                end
            end
            
            P_vectors = [P_vectors; V_P'];
            D1_vectors = [D1_vectors; V_D1'];
            D2_vectors = [D2_vectors; V_D2'];
            A_vectors = [A_vectors; A'];

  

            [bifurcation_angles]=extract_angles_lee(V_P, V_D1, V_D2, A);
%            bifurcation_angles = real(bifurcation_angles);
            for k=1:length(bifurcation_angles)
                if imag(bifurcation_angles(k))~=0
                    disp(['Imaginary angle detected at node ', num2str(hubNode), '. Angle ', num2str(i), ' is ', num2str(bifurcation_angles(k))])
                elseif isnan(bifurcation_angles(k))
                    disp(['NaN angle detected at node ', num2str(hubNode)])
                end
            end
            
            hubs = [hubs; hubNode];
            psi_1 = [psi_1; bifurcation_angles(1)];
            psi_2 = [psi_2; bifurcation_angles(2)];

            P_CPs = [P_CPs; numCPs_P];
            D1_CPs = [D1_CPs; numCPs_D1];
            D2_CPs = [D2_CPs; numCPs_D2];
            
%             % Need to shift the phi_2 angle to be between [0, 2pi). Do so by
%             % examing the x and y coordinates of the D2 vector, adjusted to the
%             % standard reference frame.
% 
%             % Standard reference frame
%             E1=[1,0,0];
%             E2=[0,1,0];
%             E3=[0,0,1];
% 
%             % Vessel trajectories
%             P=V_P/norm(V_P);
%             D1=V_D1/norm(V_D1);
%             D2=V_D2/norm(V_D2);
% 
%             % Current frame for this junction
%             e3=-P/norm(P);
%             e1=A-dot(A, e3)*e3;
%             e1=e1/norm(e1);
%             e2=cross(e3, e1);
% 
%             R=zeros([3,3]);
%             R(1,1)=dot(e1,E1);    R(1,2)=dot(e1,E2);    R(1,3)=dot(e1,E3);
%             R(2,1)=dot(e2,E1);    R(2,2)=dot(e2,E2);    R(2,3)=dot(e2,E3);
%             R(3,1)=dot(e3,E1);    R(3,2)=dot(e3,E2);    R(3,3)=dot(e3,E3);
% 
%             D2_std = R * D2;
% 
%             % Shift phi_2
%             phi_2_measure = bifurcation_angles(8);
%             if D2_std(1) < 0 
%                 phi_2_measure = phi_2_measure + pi;
%             elseif D2_std(2) < 0
%                 phi_2_measure = phi_2_measure + 2*pi;
%             end
                
            theta_1 = [theta_1; pi-bifurcation_angles(3)];
            theta_2 = [theta_2; pi-bifurcation_angles(4)];
            
            % Categorize each changepoint method
            numCPs = [numCPs_P, numCPs_D1, numCPs_D2];
            for j=1:length(numCPs)
                if numCPs(j)==0
                    n0cp = n0cp+1;
                elseif numCPs(j)==1
                    n1cp = n1cp+1;
                elseif numCPs(j)==2
                    n2cp = n2cp+1;
                    disp(['Node ', num2str(hubNode), ' has a vessel vector with 2 change points.'])
                end
            end
        end
        
    % If we are at a non-bifurcation node, we don't compute angles.
    elseif numD>=3
        msg=['Node ', num2str(hubNode), ' has ' num2str(numD) ' daughter vessels.'];
        %disp(msg);
%     elseif numD==1 && hubNode==root_node
%         msg=['Node ', num2str(hubNode), ' has 1 daughter vessel because it is the root node.'];
%         %disp(msg);
%     elseif numD==1 && hubNode~=root_node
%         msg=['Node ', num2str(hubNode), ' has 1 daughter vessel but is not the root node.'];
%         error(msg);
    end

%     if hubNode == 532
%         error('stop')
%     end
    
end

Angles{2, 1} = hubs;
Angles{2, 2} = psi_1;
Angles{2, 3} = psi_2;
Angles{2, 4} = theta_1;
Angles{2, 5} = theta_2;

CP_Groups{2,1} = hubs;
CP_Groups{2,2} = P_CPs;
CP_Groups{2,3} = D1_CPs;
CP_Groups{2,4} = D2_CPs;

CP_Counts{2,1} = n0cp;
CP_Counts{2,2} = n1cp;
CP_Counts{2,3} = n2cp;

CP_matrix{2,2} = length(find(P_CPs == 0));
CP_matrix{2,3} = length(find(P_CPs == 1));
CP_matrix{2,4} = length(find(P_CPs == 2));
CP_matrix{3,2} = length(find(D1_CPs == 0));
CP_matrix{3,3} = length(find(D1_CPs == 1));
CP_matrix{3,4} = length(find(D1_CPs == 2));
CP_matrix{4,2} = length(find(D2_CPs == 0));
CP_matrix{4,3} = length(find(D2_CPs == 1));
CP_matrix{4,4} = length(find(D2_CPs == 2));

vessel_vectors{2,1} = hubs;
vessel_vectors{2,2} = P_vectors;
vessel_vectors{2,3} = D1_vectors;
vessel_vectors{2,4} = D2_vectors;
vessel_vectors{2,5} = A_vectors;

end
