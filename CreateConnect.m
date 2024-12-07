function [newVesConn,mapIDs,TermVes] = CreateConnect(details,maxDaughters)

NumVes = length(details(2:end,1));
VesConn = zeros(NumVes,4);
for i = 2:NumVes+1
    VesConn(i-1,1) = str2num(details{i,1});
    dau = cell2mat(details{i,8});
    if maxDaughters == 2
        if length(dau) == 1 & dau(1)~=0
            VesConn(i-1,2) = dau(1);
        elseif length(dau) == 2 
            VesConn(i-1,2) = dau(1);
            VesConn(i-1,3) = dau(2);
        end
    elseif maxDaughters == 3
        if length(dau) == 1 & dau(1)~=0
            VesConn(i-1,2) = dau(1);
        elseif length(dau) == 2
            VesConn(i-1,2) = dau(1);
            VesConn(i-1,3) = dau(2);
        elseif length(dau) == 3
            VesConn(i-1,2) = dau(1);
            VesConn(i-1,3) = dau(2);
            VesConn(i-1,4) = dau(3);
       end
   end
end

newVesConn = zeros(size(VesConn));
newVesConn(:,1) = 0:length(VesConn)-1;
for i = 1:NumVes
    if VesConn(i,2) ~= 0 
       ID = find(VesConn(:,1)==VesConn(i,2));
       newVesConn(i,2) = ID-1;
    end
    if VesConn(i,3) ~= 0 
       ID = find(VesConn(:,1)==VesConn(i,3));
       newVesConn(i,3) = ID-1;
    end
    if VesConn(i,4) ~= 0 
       ID = find(VesConn(:,1)==VesConn(i,4));
       newVesConn(i,4) = ID-1;
    end
end

TermVes = find(newVesConn(:,2)==0)-1;
mapIDs  = [(0:NumVes-1)' VesConn(:,1)];