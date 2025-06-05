% This program analyzes the AIC and BIC scores of vessel vectors 
% with 0, 1, or 2 change points
% Input: xyz = m x 3 matrix of points along an edge.
% Output: V = 3 x 1 vector chosen by AIC or BIC

function [V, min_AIC, min_BIC, numCPs] = chooseVesselVector(xyz, hubNode_coordinates)

% To make this easier, make sure your hubNode is the first listed node in
% your list.
if xyz(end,:)==hubNode_coordinates
    xyz = flipud(xyz);
elseif xyz(1,:)~=hubNode_coordinates
    error('The hub node chosen does not start or end this vessel.')
end

m = size(xyz, 1); % Number of data points
AIC = [];
BIC = [];
numCPs = 0;

if m >= 3
    [~, J_LR, ~] = linReg(xyz, hubNode_coordinates, 0);
    [AIC_LR, BIC_LR] = infoCriteria(m, 1, J_LR);
    AIC = [AIC; AIC_LR];
    BIC = [BIC; BIC_LR];
    if m >= 5
        [cp, xyz1_1CP, ~, ~, J_1CP]=changePoints(xyz, hubNode_coordinates, 0);
        if ~isnan(cp(1)) 
            [AIC_1CP, BIC_1CP] = infoCriteria(m, 3, J_1CP);
            AIC = [AIC; AIC_1CP];
            BIC = [BIC; BIC_1CP];
        end
        if m >= 7
            [cp1, cp2, xyz1_2CP, ~, ~, ~, J_2CP] = twoChangePoints(xyz, hubNode_coordinates, 0);
            if ~isnan(cp1(1)) && ~isnan(cp2(1))
                [AIC_2CP, BIC_2CP] = infoCriteria(m, 5, J_2CP);
                AIC = [AIC; AIC_2CP];
                BIC = [BIC; BIC_2CP];
            end
        end
    end

%     AIC
%     BIC

    if ~isempty(AIC)
        [min_AIC, ind_AIC] = min(AIC);
        [min_BIC, ind_BIC] = min(BIC);

        if ind_AIC == ind_BIC
            numCPs = ind_AIC-1;
            %disp(['AIC and BIC scores are lowest for numCPs=' num2str(numCPs)])
        else
            numCPs = ind_BIC-1;
            %disp(['At node (' num2str(hubNode_coordinates) ') AIC selects numCPs=' num2str(ind_AIC-1) ', but BIC chooses numCPs=' num2str(numCPs)])
        end
    end

    if numCPs == 0
        [U, ~, ~] = linReg(xyz, hubNode_coordinates, 0); 
    elseif numCPs == 1
        [U, ~, ~] = linReg(xyz1_1CP, hubNode_coordinates, 0); 
        %changePoints(xyz, hubNode_coordinates, 1)
    elseif numCPs == 2
        [U, ~, ~] = linReg(xyz1_2CP, hubNode_coordinates, 0); 
        %twoChangePoints(xyz, hubNode_coordinates, 1)
    end

    V = U(:, 1);
    
elseif m==2
    V = xyz(2, :) - xyz(1, :);
    min_AIC = nan;
    min_BIC = nan;
end

end