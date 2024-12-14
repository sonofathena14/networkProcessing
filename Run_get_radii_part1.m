function Run_get_radii_part1(Name)
close all;

s1 = strcat('../',Name,'_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m_data.txt');
s2 = strcat('../',Name,'_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m.dot');
s3 = strcat('../',Name,'_DMAP.nrrd');

[arcs,nodes]=formatData(s1,s2,s3);

plotSlicerData(arcs,nodes,'b',1); % If root node is not know run plotSlicerData(arcs,nodes,'b');

prompt = "What is the root node? ";
x = input(prompt);

[path, arcsC, nodesC, correction_log]=correctionEngine(arcs,nodes,x);
[path2,arcsC2,nodesC2] = removeTooShort(nodesC,arcsC,path);
plotSlicerData(arcsC2,nodesC2,'b',3);
[path3,arcsC3,nodesC3] = blobRemover(nodesC2,arcsC2,path2);
plotSlicerData(arcsC3,nodesC3,'b',4);

[orientation,newNetwork,connectivity,arcsC2,maxDaughters,maxParrents] = directedGraphInfo(arcsC3,nodesC3, path3);

[vessel_details] = extract_geometry(newNetwork,connectivity,arcsC2,orientation);

[T] = Export_to_Excel(vessel_details,0,Name,'Aorta');

save(strcat('Output/Segmentation_',Name,'.mat'));
disp('Now run the R-code, input the correct Excel file name');
end