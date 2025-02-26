function plotgraph(Name) 
load(strcat('Output/Vessels_',Name,'.mat'),'arcsC3');
load(strcat('Output/Vessels_',Name,'.mat'),'nodesC2');

plotSlicerData_lines(arcsC3, nodesC2, 'b',1);

end
