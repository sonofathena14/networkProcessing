% Program to determine change-point along a vessel
% Input: xyz = [n x 3] matrix of coordinates of points along a vessel
%        ** n must be >= 5
% Output: cp = [1 x 3] vector of coordinates of the best change point
%         errors = [n x 6] matrix with columns:
%     [ x_cp | y_cp | z_cp | totalSSE_using_that_cp | SSE_first_segment | SSE_second_segment] 

function [cp, xyz1, xyz2, errors, minError]=changePoints(xyz, hubNode_coordinates, plot_on)

numPoints = size(xyz, 1);

if numPoints < 5
    error(['Cannot fit a changepoint to ' num2str(numPoints) ' points. Need at least 5.'])
end

% To make this easier, make sure your hubNode is the first listed node in
% your list.
if xyz(end,:)==hubNode_coordinates
    xyz = flipud(xyz);
elseif xyz(1,:)~=hubNode_coordinates
    error('The hub node chosen does not start or end this vessel.')
end

% First, obtain the SSE from using regular SVD regression
N = numPoints;
[~,SSE,~]=linReg(xyz, hubNode_coordinates, 0);
errors = [nan nan nan round(SSE/N, 16) nan nan];

% Next, place a change point at 2nd edge point from the hubNode.
% Then, move it to the next one, then the next, adding the change point
% coordinates and SSE's to the 'errors' matrix as you go.
for i=3:numPoints-2
    cp = xyz(i, :);
    xyz1 = xyz(1:i, :);
    N = size(xyz1, 1);
    xyz2 = xyz(i:end, :);
    M = size(xyz2, 1);
    [~,SSE1,~]=linReg(xyz1, xyz1(1, :), 0);
    [~,SSE2,~]=linReg(xyz2, xyz2(1, :), 0);
    %errors = [errors; cp round((SSE1+SSE2)/(N+M), 16) round(SSE1/N, 16) round(SSE2/M, 16)]
    errors = [errors; cp round((SSE1/N+SSE2/M), 16) round(SSE1/N, 16) round(SSE2/M, 16)];
end

[minError, I] = min(errors(:, 4));
cp = errors(I, 1:3);

if ~isnan(cp(1))
    [~, cp_index]=ismember(cp,xyz,'rows');
    xyz1 = xyz(1:cp_index, :);
    xyz2 = xyz(cp_index:end, :);
    if plot_on == 1
        xyz1 = xyz(1:cp_index, :);
        xyz2 = xyz(cp_index:end, :);

        hold on
        
        scatter3(xyz(:,1), xyz(:,2), xyz(:,3))
        linReg(xyz1, hubNode_coordinates, plot_on);
        linReg(xyz2, cp, plot_on);
    end
else
    xyz1 = xyz;
    xyz2 = nan;
    if plot_on == 1
        hold on
        scatter3(xyz(:,1), xyz(:,2), xyz(:,3))
        linReg(xyz1, hubNode_coordinates, plot_on);
    end
end
    
end