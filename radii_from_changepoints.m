% Code to find change points in vascular data

function [new_vessel_details] = radii_from_changepoints(Name)
    close all; 
    
    % Load vessel data 
    %sM = strcat('./Output/Segmentation_',Name,'.mat');
    sM = strcat('./Output/Vessels_',Name,'.mat');
    
    load(sM,'vessel_details')
    
    num_ves = height(vessel_details);
    
    new_vessel_details = vessel_details(:,1:8);
    new_vessel_details(:,end+1) = {[]};
    new_vessel_details(:,end+1) = {[]};
    
    new_vessel_details(1,9:10) = {'CP_COORDINATES','Coefficients'};
    %%
    for m=2:num_ves % Loop through each vessel
    
        vessel_info = cell2mat(vessel_details(m,2)); % Extract vessel info
        radii = flip(vessel_info(:,4)); % Extract radii and convert to cm
        distance= linspace(0,cell2mat(vessel_details(m,3)),length(radii)); % Create vector for distance in vessel
     
        % Call change point finder function
        [cp_index,segments,cp_num,section,slopes] = changepoint_finder(radii,distance);
       
        % Add x,y,z coordinates of change points 
        cpx_coordinates = vessel_info(cp_index,1);
        cpy_coordinates = vessel_info(cp_index,2);
        cpz_coordinates = vessel_info(cp_index,3);
    
        cp_coords = [cpx_coordinates cpy_coordinates cpz_coordinates];
    
        new_vessel_details(m,9) = {cp_coords};
    
        % Slopes of the lines between the change points
        new_vessel_details(m,10) = {slopes'};
    end
    save(strcat('./Input/Changepoints_',Name,'.mat'));
end
    