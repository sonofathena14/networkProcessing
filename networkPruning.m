function networkPruning(Name,avg_control_branches)
% clear all
% Name = 'm1p4_053007';
% avg_control_branches = 850;

oldData = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'Data');
oldData = oldData.('Data');
oldSegments = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'segments');
oldSegments = oldSegments.('segments');
arcsC4 = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'arcsC3');
arcsC4 = arcsC4.('arcsC3');
nodesC3 = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'nodesC2');
nodesC3 = nodesC3.('nodesC2');
sf = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'sf');
sf = sf.('sf');

mapIDs = oldData.mapIDs;
TaperVes = [];
ploton = 0;
[vessels_to_prune, radThreshold]=radiusPruning(oldData, avg_control_branches);

%% remove pruned vessels
for i=1:length(vessels_to_prune)
    oldID = mapIDs(find(mapIDs(:,1)==vessels_to_prune(i)),2);
    from = oldSegments(oldID,2);
    to = oldSegments(oldID,3);
    [arcsC4,nodesC3] = edgePrune(arcsC4,nodesC3,from,to);
end


%% Reextract radii, angles and then save updated network
%calls necessary functions, check specific files for details
[orientation, newNetwork,connectivity,arcsC5,maxDaughters,~] = directedGraphInfoPrune(arcsC4,nodesC3);

[vessel_details] = extract_geometry(newNetwork,connectivity,arcsC5,orientation);

[new_vessel_details] = radii_from_changepoints(vessel_details);
    
changepoint_location=get_changepoint_locations(new_vessel_details);

vessel_radii=get_radii_mice(new_vessel_details, changepoint_location,TaperVes,ploton);

[Angles, CP_Groups, CP_Counts, CP_matrix, vessel_vectors] = extractNewNetworkAngles_MAIN(arcsC5, nodesC3, newNetwork, connectivity);

segments=zeros(length(arcsC5),3); %creates a list of ves id with from and to nodes
for i=1:length(arcsC5)
    segments(i,1)=i;
    segments(i,2)=arcsC5{1,i}(1,1);
    segments(i,3)=arcsC5{1,i}(1,2);
end

volumes = edgeVolume(vessel_details);

Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);
save(strcat('Pruned/Pruned_Network_',Name,'.mat'), 'arcsC5', 'nodesC3', 'Data','vessel_details', 'TaperVes','ploton','maxDaughters','segments','volumes','sf','Angles','avg_control_branches')
    % ^ saves the necessary data for statistic extraction and to run/work
    % on the changepoint algorithm

end