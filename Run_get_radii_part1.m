function Run_get_radii_part1(Name)
close all;
%Name = 'control_3';
s1 = strcat('../',Name,'_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m_data.txt');
s2 = strcat('../',Name,'_SKEL_adj_dmax_isthmus1_p0_REDUCED_sp_c_m.dot');
s3 = strcat('../',Name,'_DMAP.nrrd');

[arcs,nodes] = formatData(s1,s2,s3);
[arcs]       = OrientArcs(arcs,nodes);

plotSlicerData(arcs,nodes,'b',1); % If root node is not know run plotSlicerData(arcs,nodes,'b');
prompt = "What is the root node? ";
x = input(prompt);

[path, arcsC, nodesC, correction_log]=correctionEngine(arcs,nodes,x);
[path2,arcsC2,nodesC2] = removeTooShort(nodesC,arcsC,path);
plotSlicerData(arcsC2,nodesC2,'b',2);
i = 3;
ans1 = input('Do you need to remove any blobs?','s');
if ans1 == 'Y'
    STILL_GOING = 'Y';
    while STILL_GOING =='Y'
        [path2,arcsC2,nodesC2] = blobRemover(nodesC2,arcsC2,path2);
        plotSlicerData(arcsC2,nodesC2,'b',i);
        STILL_GOING=input('Do you need to remove another blob?','s');
        i = i+1;
    end
end

plotSlicerData_lines(arcsC2,nodesC2,'b',i);
i = i+1;
notDone = 'Y';
while notDone == 'Y'
    ans2=input('Do you need to remove any edges?','s');
    if ans2 == 'Y'
        STILL_GOING = 'Y';
        while STILL_GOING =='Y'
            [arcsC2,nodesC2,path2]=edgeRemoval(arcsC2,nodesC2,path2);
            plotSlicerData(arcsC2,nodesC2,'b',i);
            STILL_GOING=input('Do you need to remove another edge?','s');
            i = i+1;
        end
    end
    
    ans3=input('Do you need to remove any nodes?','s');
    if ans3 == 'Y'
        STILL_GOING = 'Y';
        while STILL_GOING =='Y'
            [arcsC2,nodesC2,path2]=nodeRemoval(arcsC2,nodesC2,path2);
            plotSlicerData(arcsC2,nodesC2,'b',i);
            STILL_GOING=input('Do you need to remove another node?','s');
            i = i+1;
        end
    end
    notDone=input('Do you need to continue (Y/N)?','s');
end


[orientation,newNetwork,connectivity,arcsC3,maxDaughters,maxParrents] = directedGraphInfo(arcsC2,nodesC2, path2);

[vessel_details] = extract_geometry(newNetwork,connectivity,arcsC3,orientation);
save(strcat('Output/Vessels_',Name,'.mat'));

% [T] = Export_to_Excel(vessel_details,0,Name,'Aorta');
% save(strcat('Output/Segmentation_',Name,'.mat'));

disp('Now run next code, ensuring the name is correct');

end
