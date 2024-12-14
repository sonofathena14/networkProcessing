function [new_vessel_details,Data,newNetwork,DataOut] = Run_get_radii_part2(Name,TaperVes,Scale,ploton)
    sM = strcat('./Output/Segmentation_',Name,'.mat');
    sE = strcat('./Input/IMPORT_',Name,'.xlsx');
    sO = strcat('./FluidsInput/FluidInput_',Name,'.mat');

    load(sM,'maxDaughters','vessel_details','newNetwork','arcsC2','nodesC');
    
    [new_vessel_details] = Import_from_excel(sE,vessel_details);
    
    changepoint_location=get_changepoint_locations(new_vessel_details);

    disp('ttt')
    vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

    disp('done with get_radii')
    Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);

    % [DataOut] = get_radii_aorta(new_vessel_details, Scale, Data);  % Only call if you want to correct radii extraction for tapering aorta.
    DataOut = Data;
    save(sO,'Data','DataOut');
end