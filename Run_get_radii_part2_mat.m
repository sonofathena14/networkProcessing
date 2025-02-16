% load
% vessel_details

function [new_vessel_details,Data,newNetwork] = Run_get_radii_part2_mat(Name,TaperVes,Scale,ploton)

% TaperVes = [];
% ploton = 1;
% Scale  = 0.742;

    sM = strcat('./Output/Segmentation_',Name,'.mat');
    sE = strcat('./Input/Changepoints_',Name,'.mat'); %replace with correct matlab import 
    sO = strcat('./FluidsInput/FluidInput_',Name,'.mat');
    load(sM,'maxDaughters','newNetwork');
    load(sE,'new_vessel_details')
    
    changepoint_location=get_changepoint_locations(new_vessel_details);
    vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

    Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);
    save(sO,'Data');
end

    
    