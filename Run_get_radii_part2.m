function [new_vessel_details,Data,newNetwork] = Run_get_radii_part2(Name,TaperVes,Scale,ploton)

    clear all;
    close all;

    %Name = 'control_3'; % name of network
    % TaperVes= []; % index of vessels to be tapered put in "XX"
    % Scale = 1;   % scaling from voxels to mm from CT image
    % ploton = 1;  % plot figures

    sM = strcat('./Output/Segmentation_',Name,'.mat');
    sE = strcat('./Input/IMPORT_',Name,'.xlsx'); 
    sO = strcat('./FluidsInput/FluidInput_',Name,'.mat');

    load(sM,'maxDaughters','vessel_details','newNetwork','arcsC3','nodesC2');
    
    [new_vessel_details] = Import_from_excel(sE,vessel_details);
    

    changepoint_location=get_changepoint_locations(new_vessel_details);
    save(strcat('FluidsInput/Changepoint_',Name,'.mat'))

    vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

    disp('done with get_radii')
    Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);

    %[DataOut] = get_radii_aorta(new_vessel_details, Scale, Data);  % Only call if you want to correct radii extraction for tapering aorta.

    save(sO,'Data'); %only saves the necessary input for the fluids code
    save(strcat('Output/Vessels_',Name,'.mat')); %saves entire workspace
    save(strcat('Networks/Network_Vessels_',Name,'.mat'), 'arcsC3', 'nodesC2', 'Data','vessel_details', 'TaperVes','ploton','Scale','maxDaughters')
        % ^ saves the necessary data for statistic extraction and to run/work
        % on the changepoint algorithm
end