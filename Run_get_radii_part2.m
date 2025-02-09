function [new_vessel_details,Data,newNetwork] = Run_get_radii_part2(Name,TaperVes,Scale,ploton)
    % Name = 'PVB046_arteries'; %name of network
    % TaperVes=[]; index of vessels to be tapered
    % Scale = 1; scaling from voxels to mm from CT image
    % ploton = 1; 1 means plot figures
    
    sM = strcat('./Output/Segmentation_',Name,'.mat');
    sE = strcat('./Input/IMPORT_',Name,'.xlsx');
    sO = strcat('./FluidsInput/FluidInput_',Name,'.mat');

    load(sM,'maxDaughters','vessel_details','newNetwork','arcsC2','nodesC');
    
    [new_vessel_details] = Import_from_excel(sE,vessel_details);
    
    changepoint_location=get_changepoint_locations(new_vessel_details);
    save(strcat('FluidsInput/Changepoint_',Name,'.mat'))
    
    vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

    Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);
    save(sO,'Data','new_vessel_details');
end