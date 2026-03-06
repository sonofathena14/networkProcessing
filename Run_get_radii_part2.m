function [new_vessel_details,Data,newNetwork] = Run_get_radii_part2(Name,TaperVes,Scale,ploton)

    clear all;
    close all;

    Name = 'AortaTest'; % name of network
    TaperVes= []; % index of vessels to be tapered put in "XX"
    Scale = 1.041;   % scaling from voxels to mm from CT image
    ploton = 1;  % plot figures

    fluids_data_table_filename = Name + "_Fluids_Data_Table.mat";

    if isfile(fluids_data_table_filename)
    load(fluids_data_table_filename, 'T_fluids');
    end

    sM = strcat('./Output/Segmentation_',Name,'.mat');
    sE = strcat('./Input/IMPORT_',Name,'.xlsx'); 
    sO = strcat('./FluidsInput/FluidInput_',Name,'.mat');

    recycle on;
    if isfile(sO)
        delete(sO);
    end


    load(sM,'maxDaughters','vessel_details','newNetwork','arcsC3','nodesC2','connectivity');
    
    [new_vessel_details] = Import_from_excel(sE,vessel_details);
    

    changepoint_location=get_changepoint_locations(new_vessel_details);
    %save(strcat('FluidsInput/Changepoint_',Name,'.mat'))

    vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

    [Angles, CP_Groups, CP_Counts, CP_matrix, vessel_vectors] = extractNewNetworkAngles_MAIN(arcsC3, nodesC2, newNetwork, connectivity);

    segments=zeros(length(arcsC3),3); %creates a list of ves id with from and to nodes
    for i=1:length(arcsC3)
        segments(i,1)=i;
        segments(i,2)=arcsC3{1,i}(1,1);
        segments(i,3)=arcsC3{1,i}(1,2);
    end

    
    volumes = edgeVolume(vessel_details);

    plotSlicerData(arcsC3,nodesC2,'b',1);

    disp('done with get_radii')
    Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters,segments);

    assignin('base', ['vessel_rad' ...
        'ii'], vessel_radii);
    [DataOut] = get_radii_aorta_old(new_vessel_details, Scale, Data);  % Only call if you want to correct radii extraction for tapering aorta.

    % Pad or truncate mapIDs to fit 20 rows
numMapRows = size(DataOut.mapIDs, 1);
if numMapRows <= 20
    T_fluids.Map_ID_From(1:numMapRows) = DataOut.mapIDs(:,1);
    T_fluids.Map_ID_To(1:numMapRows) = DataOut.mapIDs(:,2);
else
    warning('mapIDs has %d rows but table has 20. Using first 20 rows.', numMapRows);
    T_fluids.Map_ID_From = DataOut.mapIDs(1:20,1);
    T_fluids.Map_ID_To = DataOut.mapIDs(1:20,2);
end


save(fluids_data_table_filename, 'T_fluids');

    save(sO,'Data','DataOut'); %only saves the necessary input for the fluids code
    save(strcat('Output/Vessels_',Name,'.mat')); %saves entire workspace
    save(strcat('Networks/Network_Vessels_',Name,'.mat'), 'arcsC3', 'nodesC2', 'Data','vessel_details', 'TaperVes','ploton','Scale','maxDaughters','Angles','volumes');
        % ^ saves the necessary data for statistic extraction and to run/work
        % on the changepoint algorithm
end