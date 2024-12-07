function [T] = Export_to_Excel(details,Extraction,FN,Name)

%% automate for entire dataset
% extract the correct information from each vessel (length, radius, xyz
% coordinates, parent & daughter(s)). 
% If centerlines are extracted by VMTK set Extraction = 1 else set Extraction = 0.
% This is needed since we have to flip the radius
% and xyz components because VMTK works backwards from the outlets to the
% inlet(s).

excelname = strcat('EXPORT_',FN,'.xlsx');
     
for i = 2:length(details(:,1))
    if (mod(i,100) == 0) disp([i]); end
    if Extraction == 1 % Flip radii
        Radius       = flip(details{i,2}(2:end,4)./10); %converted from mm to cm
        Length       = [linspace(0,details{i,3},length(Radius))]'; %converted from mm to cm
        x            = flip(details{i,2}(2:end,1));
        y            = flip(details{i,2}(2:end,2));
        z            = flip(details{i,2}(2:end,3));
    else
        Radius       = details{i,2}(2:end,4)./10; %converted from mm to cm
        Length       = [linspace(0,details{i,3},length(Radius))]'; %converted from mm to cm Note -  % Length is not necessarily equidistant between points (there should be a length segment between each xyz coordinate)
        x            = details{i,2}(2:end,1);
        y            = details{i,2}(2:end,2);
        z            = details{i,2}(2:end,3);
    end
    VesID            = details{i,1};
    
    Sheet_Name = strcat(Name,'_',VesID); %Takes the name of the vessel and makes it the corresponding sheet name in excel
    T  = table(Radius, Length, x, y, z); %table with all numerical values
   
    % Export to excel where each vessel has it's own spreadsheet. The columns
    % are in the following order: Radius, Length, xyz, Daughter 1, 2, etc.,
    % Parent. In R two additional columns will be added that are the location
    % of the changepoints and the slopes of the piecewise functions,
    % respectively. The excel file is written to the current working directory
    % unless otherwise specified.
    
    writetable(T,excelname,'Sheet',Sheet_Name);
end 
end