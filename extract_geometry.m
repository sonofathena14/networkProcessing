function [vessel_details] = extract_geometry(newNetwork, connectivity, arcsC, orientation)

rootIDrow    = find(connectivity(:, 3)==0);
rootID = connectivity(rootIDrow,1);
rootVesID = find(newNetwork(:,2) == rootID);

NumVes    = length(newNetwork);
AllVesID  = find(newNetwork (:,1));

nodeOUT= newNetwork(:,3);
nodeIN = newNetwork(:,2);
TermVesID = newNetwork(~ismember(nodeOUT, nodeIN), 1);
IntVesID  = setdiff(AllVesID,[TermVesID; rootVesID]);

% Create vessel details, structure is given below
vessel_details = cell(NumVes, 9);
vessel_details{1, 1} = 'VES ID';
vessel_details{1, 2} = 'CL FILE';
vessel_details{1, 3} = 'LENGTH [mm]';
vessel_details{1, 4} = 'MEAN RADII [mm]';
vessel_details{1, 5} = 'MEDIAN RADII';
vessel_details{1, 6} = 'R_STD';
vessel_details{1, 7} = 'PARENT';
vessel_details{1, 8} = 'DAUGHTERS';
vessel_details{1, 9} = 'DAUGHTERS ID'; %%%hell this new

% Root vessel
vessel_details {2,1} = num2str(rootVesID);
arcsC2 = arcsC{1,rootVesID};
Nrow = size(arcsC{1,rootVesID},1);
if orientation(rootVesID) == -1
    arcsC2(1,1)    = arcsC{1,rootVesID}(1,2);
    arcsC2(1,2)    = arcsC{1,rootVesID}(1,1);
    arcsC2(2:Nrow,1:4)= flipud(arcsC{1,rootVesID}(2:Nrow,1:4));
end
vessel_details {2,2} = arcsC2;
vessel_details {2,3} = newNetwork (rootVesID,5);
vessel_details {2,4} = mean(arcsC{1, rootVesID}(2:end,4));
vessel_details {2,5} = median(arcsC{1, rootVesID}(2:end,4));
vessel_details {2,6} = std(arcsC{1, rootVesID}(2:end,4));
vessel_details {2,7} = 0;
DauVesID = find(newNetwork(:,2)==newNetwork(rootVesID,3));
vessel_details {2,8} = {DauVesID};

if length(DauVesID) == 2
    vessel_details{2,9} = strcat(num2str(DauVesID(1)), ',', num2str(DauVesID(2)));
elseif length(DauVesID) == 3
    vessel_details{2,9} = strcat(num2str(DauVesID(1)), ',', num2str(DauVesID(2)), ',', num2str(DauVesID(3)));
elseif length(DauVesID) == 4
    vessel_details{2,9} = strcat(num2str(DauVesID(1)), ',', num2str(DauVesID(2)), ',', num2str(DauVesID(3)), ',', num2str(DauVesID(4)));
end

NtermVes = length(TermVesID);
for i = 1:NtermVes % Terminal vessels
    vessel_details {i+2,1} = num2str(TermVesID(i));
    arcsC2 = arcsC{1,TermVesID(i)};
    Nrow = size(arcsC{1,TermVesID(i)},1); 
    if orientation(TermVesID(i)) == -1
        arcsC2(1,1) = arcsC{1,TermVesID(i)}(1,2);
        arcsC2(1,2) = arcsC{1,TermVesID(i)}(1,1);
        arcsC2(2:Nrow,1:4) = flipud(arcsC{1,TermVesID(i)}(2:Nrow,1:4));
    end

    vessel_details {i+2,2} = arcsC2;
    vessel_details {i+2,3} = newNetwork (TermVesID(i),5);
    vessel_details {i+2,4} = mean(arcsC{1,TermVesID(i)}(2:end,4));
    vessel_details {i+2,5} = median(arcsC{1,TermVesID(i)}(2:end,4));
    vessel_details {i+2,6} = std(arcsC{1,TermVesID(i)}(2:end,4));
    inNode = newNetwork(TermVesID(i),2);
    parVesID = find(newNetwork(:, 3) == inNode);
    vessel_details {i+2,7} = parVesID;
    vessel_details {i+2,8} = {0};
    vessel_details{i+2,9} = strcat(num2str(0), ',', num2str(0),',', num2str(0));
end

for i = 1:length(IntVesID) % Interior vessels
    vessel_details {i+2+NtermVes,1} = num2str(IntVesID(i));
    arcsC2 = arcsC{1,IntVesID(i)};
    Nrow = size(arcsC{1,IntVesID(i)},1);
    if orientation(IntVesID(i)) == -1
        arcsC2(1,1) = arcsC{1,IntVesID(i)}(1,2);
        arcsC2(1,2) = arcsC{1,IntVesID(i)}(1,1);
        arcsC2(2:Nrow,1:4) = flipud(arcsC{1,IntVesID(i)}(2:Nrow,1:4));
    end
    vessel_details {i+2+NtermVes,2} = arcsC2;
    vessel_details {i+2+NtermVes,3} = newNetwork (IntVesID(i),5);
    vessel_details {i+2+NtermVes,4} = mean(arcsC{1,IntVesID(i)}(2:end,4));
    vessel_details {i+2+NtermVes,5} = median(arcsC{1,IntVesID(i)}(2:end,4));
    vessel_details {i+2+NtermVes,6} = std(arcsC{1,IntVesID(i)}(2:end,4));
    inNode = newNetwork(IntVesID(i), 2);
    parVesID = find(newNetwork(:, 3) == inNode);
    vessel_details {i+2+NtermVes,7} = parVesID ;

    daughters = connectivity(connectivity(:, 8) == IntVesID(i), 4:7);
    NonzeroDauVesID = daughters(daughters ~= 0);

    if ~isempty(NonzeroDauVesID)
        DauVesID = NonzeroDauVesID;
        vessel_details {i+2+NtermVes,8} = {DauVesID};

        if length(DauVesID) == 2
            vessel_details{i+2+NtermVes,9} = strcat(num2str(DauVesID(1)), ',', num2str(DauVesID(2)));
        elseif length(DauVesID) == 3
            vessel_details{i+2+NtermVes,9} = strcat(num2str(DauVesID(1)), ',', num2str(DauVesID(2)), ',', num2str(DauVesID(3)));
        elseif length(DauVesID) == 4
            vessel_details{i+2+NtermVes,9} = strcat(num2str(DauVesID(1)), ',', num2str(DauVesID(2)), ',', num2str(DauVesID(3)), ',', num2str(DauVesID(4)));
        end
    end
end
end 