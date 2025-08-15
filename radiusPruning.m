% This function performs radius cropping on a hypoxic arterial network.

% Pairs of terminal vessels are removed if both their radii fall below
% a certain threshold. When a pair of vessels is removes, its parent 
% becomes a terminal vessel and the algorithm continues The threshold 
% adjusts until the pruned network has a total number of branches that is 
% closest to the average number of branches in control mice of that 
% pressure (greater than or less than the average branch count - take 
% whichever is closest).

% Inputs: Data (structure), avg_control_branches (number)

% Outputs: vessels_to_prune (vector), radThreshold (number)

function [vessels_to_prune, radThreshold]=radiusPruning(Data, avg_control_branches)

%% Get minimum radius and starting threshold
radii=Data.rin;
minRad=min(radii);
radThreshold=ceil(minRad);
% Define how much you will raise radThreshold by for each iteration.
increment = 1; 

%% Make a list of vessels to remove from the tree
TV = Data.TermVes;
Conn = Data.connectivity;
Total_Vessels = length(Conn);
vessels_to_prune = [];
vessels_remaining = Total_Vessels;
while vessels_remaining >= avg_control_branches
    % Save the old list so that you can compare the number of branches
    % closest to and greater than / less than the average control branches
    vessels_to_prune_OLD = vessels_to_prune;
    vessels_remaining_OLD = vessels_remaining;

    % Search through TermVes for pairs to remove.
    i=1;
    while i<=length(TV)
        % Get TermVes ID
        TV_ID=TV(i);
        if TV_ID>=0
            % Find the row R in connectivity where TV_ID is a daughter.
            [R, ~] = find(Conn(:,2:end)==TV_ID);

            % Find parent and all TermVes IDs that are daughters here.
            daughters = Conn(R,2:end);
            daughters = daughters(daughters~=0);
            
            % Check if ALL daughters here are also terminal
            Lia = ismember(daughters',TV);
            if all(Lia)
                % All daughters are terminal here.
                % Check if radii of daughters is < radThreshold.
                num_daughters_less_than_radThreshold = 0;
                for j=1:length(daughters)
                    DR=radii(daughters(j)+1);
                    if DR < radThreshold
                        num_daughters_less_than_radThreshold = num_daughters_less_than_radThreshold+1;
                    else
                        break
                    end
                end

                % If num_daughters_less_than_radThreshold = length(daughters),
                % add the daughters to the vessels_to_prune list.
                if num_daughters_less_than_radThreshold == length(daughters)
                    % Add daughters to vessels_to_prune list.
                    vessels_to_prune = [vessels_to_prune; daughters'];
                    % Remove these from TermVes list. This is done by
                    % replacing their IDs in TermVes with -1.
                    for j=1:length(daughters)
                        % Replace these IDs in TV with -1 to indicate these
                        % terminal vessels are removed.
                        TV(find(TV(:)==daughters(j)))=-1;
                    end
                    
                    % Add the parent of these terminal vessels to TV, 
                    % since it is now terminal.
                    TV = [TV; Conn(R, 1)];
                end
            end
        end

        % Next entry
        i=i+1;
    end

    % Recalculate vessels_remaining
    vessels_remaining = Total_Vessels - length(vessels_to_prune);

    % Calculate new radThreshold - raising by increment.
    radThreshold=radThreshold+increment;
end

%% Decide whether to prune above/below avg_control_branches

% avg_control_branches - pruned_tree_branches if pruned is BIGGER
difference1 = vessels_remaining_OLD - avg_control_branches;

% avg_control_branches - pruned_tree_branches if pruned is SMALLER
difference2 =  - avg_control_branches - vessels_remaining;

if difference1 <= difference2
    % The pruned branch count was closer in the BIGGER pruned tree than it
    % was in the SMALLER pruned tree. So go back 1 iteration to that tree.
    vessels_to_prune = vessels_to_prune_OLD;
    vessels_remaining=vessels_remaining_OLD;
    radThreshold=radThreshold-increment;
end

disp(['Pruned tree contains ', num2str(vessels_remaining), ' vessels.'])
disp(['That is ', num2str(vessels_remaining-avg_control_branches), ...
    ' above the average control branch count of ', num2str(avg_control_branches),'.'])
disp(['The radius threshold for this pruned tree is ', num2str(radThreshold), ' microns.'])

end