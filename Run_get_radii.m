function Run_get_radii(Name,TaperVes,Scale,ploton)
close all;
% Name = 'm2p4_053007';
% TaperVes = [];
% Scale = 1;
% ploton=0;
s1 = strcat('../',Name,'_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m_data.txt');
s2 = strcat('../',Name,'_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m.dot');
s3 = strcat('../',Name,'_DMAP.nrrd');
sO = strcat('./FluidsInput/FluidInput_',Name,'.mat');

[arcs,nodes] = formatData(s1,s2,s3);
[arcs]       = OrientArcs(arcs,nodes);

plotSlicerData(arcs,nodes,'b',1); %If root node is not know run plotSlicerData(arcs,nodes,'b',1);
prompt = "What is the root node? ";
x = input(prompt);

[path, arcsC, nodesC, correction_log]=correctionEngine(arcs,nodes,x);
plotSlicerData(arcsC,nodesC,'b',2);
[path2,arcsC2,nodesC2] = removeTooShort_norm(nodesC,arcsC,path,x); %removes any terminal vessel shorter than 5 voxels
plotSlicerData(arcsC2,nodesC2,'b',3);
i = 4;
ans1 = input('Do you need to remove any blobs?','s');
if ans1 == 'Y'
    STILL_GOING = 'Y';
    while STILL_GOING =='Y'
        [path2,arcsC2,nodesC2] = blobRemover(nodesC2,arcsC2,path2,x);
        plotSlicerData(arcsC2,nodesC2,'b',i);
        STILL_GOING=input('Do you need to remove another blob?','s');
        i = i+1;
    end
end

%plotSlicerData_lines(arcsC2,nodesC2,'b',i);
plotSlicerData(arcsC2,nodesC2,'b',i);
i = i+1;
trians = input('Do you need to generate trifurcations?','s');
if trians == 'Y'
  STILL_GOING = 'Y'; 
  while STILL_GOING =='Y'
     [path2,arcsC2,nodesC2] = generateTrifurcation(nodesC2,arcsC2,path2,x);
     plotSlicerData(arcsC2,nodesC2,'b',i);
     STILL_GOING = input('Do you need to generate another trifurcation?','s');
     i = i+1;
  end
end

notDone = 'Y';
while notDone == 'Y'
    ans2=input('Do you need to remove any edges?','s');
    if ans2 == 'Y'
        STILL_GOING = 'Y';
        while STILL_GOING =='Y'
            [arcsC2,nodesC2,path2]=edgeRemoval(arcsC2,nodesC2,path2,x);
            plotSlicerData(arcsC2,nodesC2,'b',i);
            STILL_GOING=input('Do you need to remove another edge?','s');
            i = i+1;
        end
    end
    notDone=input('Do you need to continue (Y/N)?','s');
end

%calls necessary functions, check specific files for details
[orientation,newNetwork,connectivity,arcsC3,maxDaughters,maxParrents] = directedGraphInfo(arcsC2,nodesC2, path2);

[vessel_details] = extract_geometry(newNetwork,connectivity,arcsC3,orientation);

[new_vessel_details] = radii_from_changepoints(vessel_details);
    
changepoint_location=get_changepoint_locations(new_vessel_details);

vessel_radii=get_radii(new_vessel_details, changepoint_location,TaperVes,Scale,ploton);

[Angles, CP_Groups, CP_Counts, CP_matrix, vessel_vectors] = extractNewNetworkAngles_MAIN(arcsC3, nodesC2, newNetwork, connectivity);

segments=zeros(length(arcsC3),3); %creates a list of ves id with from and to nodes
for i=1:length(arcsC3)
    segments(i,1)=i;
    segments(i,2)=arcsC3{1,i}(1,1);
    segments(i,3)=arcsC3{1,i}(1,2);
end

volumes = edgeVolume(vessel_details);

Data = CreateFluidsCodeInput(vessel_radii,new_vessel_details,maxDaughters);
Data = get_radii_aorta(new_vessel_details,Data,Scale);
save(sO,'Data'); %only saves the necessary input for the fluids code
%save(strcat('Output/Vessels_',Name,'.mat')); %saves entire workspace
save(strcat('Networks/Network_Vessels_',Name,'.mat'), 'arcsC3', 'nodesC2', 'Data','vessel_details', 'TaperVes','ploton','Scale','maxDaughters','segments','volumes','Angles')
    % ^ saves the necessary data for statistic extraction and to run/work
    % on the changepoint algorithm

end
