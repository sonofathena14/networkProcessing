% This function will take a vessel dataset and as input, and will return a 
% cell that contains the name of the vessel, a matrix containing the values
% of the radius and length at each changepoint, and the slopes between each
% changepoint

% IN ORDER TO RUN THIS FUNCTION SUCCESSFULLY, YOU MUST HAVE PREVIOUSLY RUN
% ImportBack.m 

function changepoint_locations = get_changepoint_locations(details)

changepoint_locations = cell(height(details), 3);

changepoint_locations{1, 1} = "VESSEL";
changepoint_locations{1, 2} = "RAD + LEN";
changepoint_locations{1, 3} = "SLOPES";

for i = 2:height(details)
    changepoint_locations{i, 1} = details{i, 1};
    if ~isempty(details{i,9})
       num_cpts = length(details{i,9}(:,1));
    else
       num_cpts = 0;
    end
    if num_cpts > 0
      cpts = details{i,9};
      cpt_radii   = zeros(num_cpts, 1);
      cpt_lengths = zeros(num_cpts, 1);
    else
      cpts = [];
      cpt_radii = [];
      cpt_lengths = [];
    end
    
    % If vessel has cpts, find the index of the x-value of the cpt in 
    % vessel_details_table, and get the radius value and length value at
    % that location
    if ~isempty(cpts)
        %disp(strcat('I = ',num2str(i)))
        for j = 1:num_cpts
            cpt_x_value  = cpts(j, 1);
            cpt_y_value  = cpts(j, 2);
            cpt_z_value  = cpts(j, 3);
            radii_values = details{i, 2};
            ID = find(radii_values(:,1)==cpt_x_value & radii_values(:,2)==cpt_y_value & radii_values(:,3)==cpt_z_value);
            %disp(j)
            cpt_radii  (j, 1) = radii_values(ID,4);
            cpt_lengths(j, 1) = radii_values(ID,5);
        end
    end
   
    % If vessel has cpts, add them to changepoint_locations, else add empty
    % vectors
    if ~isempty(cpts)
         changepoint_locations{i, 2} = [cpt_radii, cpt_lengths];
    else
         changepoint_locations{i, 2}  = [];
    end

    % Get slopes between cpts, run regression from beginning to first cpt
    % and last cpt to end to get slope
    changepoint_locations{i,3} = details{i,10};
end

end