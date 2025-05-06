
Name = 'm3p4_060407';

load(strcat('Networks/Network_Vessels_',Name,'.mat'), 'vessel_details', 'TaperVes','Scale', 'ploton','maxDaughters');

[new_vessel_details] = radii_from_changepoints(vessel_details);
    
changepoint_location=get_changepoint_locations(new_vessel_details);

vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);

