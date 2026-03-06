function Run_get_radii_part1(Name)
close all;

geoPath = '/Users/sareyesr/SGext/build/Geometries/';

Name = 'AortaTest';

folder_name = Name + "GraphFigures";

if ~exist(folder_name,'dir')
    mkdir(folder_name);
end

fluids_data_table_filename = Name + "_Fluids_Data_Table.mat";

s1 = fullfile(geoPath, [Name, '_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m_data.txt']);
s2 = fullfile(geoPath, [Name, '_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m.dot']);
s3 = fullfile(geoPath, [Name, '_DMAP.nrrd']);

[arcs,nodes] = formatData(s1,s2,s3);
[arcs]       = OrientArcs(arcs,nodes);

plotSlicerData(arcs,nodes,'b',1);
savefig(gcf, fullfile(folder_name,'Original_Geometry.fig'));

prompt = "What is the root node? ";
x = input(prompt);

[path, arcsC, nodesC, correction_log]=correctionEngine(arcs,nodes,x);
plotSlicerData(arcsC,nodesC,'b',2);
[path2,arcsC2,nodesC2] = removeTooShort_norm(nodesC,arcsC,path,x);
plotSlicerData(arcsC2,nodesC2,'b',3);
savefig(gcf, fullfile(folder_name,'Identify_Root_Geometry.fig'));


i = 4;
m1 = 1;
notDone = 1;
while notDone == 1
    ans1 = input('Do you need to remove any blobs? Y(1) / N(0): ');
    if ans1 == 1
        STILL_GOING = 1;
        while STILL_GOING == 1
            [path2,arcsC2,nodesC2] = blobRemover(nodesC2,arcsC2,path2,x);
            plotSlicerData(arcsC2,nodesC2,'b',i);
            savefig(gcf, fullfile(folder_name,"Blob_Removal_Geometry_" + m1 + ".fig"));
            STILL_GOING = input('Do you need to remove another blob? Y(1) / N(0): ');
            i = i+1;
            m1 = m1 + 1;
        end
    end
    
    plotSlicerData(arcsC2,nodesC2,'b',i);
    
    i = i+1;
    m3 = 1;
    trians = input('Do you need to generate trifurcations? Y(1) / N(0): ');
    if trians == 1
      STILL_GOING = 1; 
      while STILL_GOING == 1
         [path2,arcsC2,nodesC2] = generateTrifurcation(nodesC2,arcsC2,path2,x);
         plotSlicerData(arcsC2,nodesC2,'b',i);
         savefig(gcf, fullfile(folder_name,"Trifurcation_Generation_" + m3 + ".fig"));
         STILL_GOING = input('Do you need to generate another trifurcation? Y(1) / N(0): ');
         i = i+1;
         m3 = m3 + 1;
      end
    end
    m2 = 1;
    ans2 = input('Do you need to remove any edges? Y(1) / N(0): ');
    if ans2 == 1
        STILL_GOING = 1;
        while STILL_GOING == 1
            [arcsC2,nodesC2,path2]=edgeRemoval(arcsC2,nodesC2,path2,x);
            plotSlicerData(arcsC2,nodesC2,'b',i);
            savefig(gcf, fullfile(folder_name,"Edge_Removal_Geometry_" + m2 + ".fig"));
            STILL_GOING = input('Do you need to remove another edge? Y(1) / N(0): ');
            i = i+1;
            m2 = m2 + 1;
        end
    end
    notDone = input('Do you need to continue? Y(1) / N(0): ');
end


[orientation,newNetwork,connectivity,arcsC3,maxDaughters,maxParents] = directedGraphInfo(arcsC2,nodesC2, path2);

[vessel_details] = extract_geometry(newNetwork,connectivity,arcsC3,orientation);

[T] = Export_to_Excel(vessel_details,0,Name,'Aorta');
save(strcat('Output/Segmentation_',Name,'.mat'));

if isfile(fluids_data_table_filename)
    load(fluids_data_table_filename, 'T_fluids');
else
    T_fluids = table( ...
        'Size', [20 5], ...
        'VariableTypes', repmat({'double'}, 1, 5), ...
        'VariableNames', {'VID', 'NewNetwork_From', 'NewNetwork_To', 'Map_ID_From', 'Map_ID_To'});
    T_fluids(:,:) = {NaN};
end

numDataRows = size(newNetwork, 1);
if numDataRows <= 20
    T_fluids.VID(1:numDataRows) = newNetwork(:,1);
    T_fluids.NewNetwork_From(1:numDataRows) = newNetwork(:,2);
    T_fluids.NewNetwork_To(1:numDataRows) = newNetwork(:,3);
else
    warning('newNetwork has %d rows but table has 20. Using first 20 rows.', numDataRows);
    T_fluids.VID = newNetwork(1:20,1);
    T_fluids.NewNetwork_From = newNetwork(1:20,2);
    T_fluids.NewNetwork_To = newNetwork(1:20,3);
end

save(fluids_data_table_filename, 'T_fluids');

save InputChangepoint.mat vessel_details;

save(strcat('Output/Segmentation_',Name,'.mat'));
disp('Now run the R-code, input the correct Excel file name');

end