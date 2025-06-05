% Program to determine change-point along a vessel
% Input: xyz = [n x 3] matrix of coordinates of points along a vessel
% Output: cp = [1 x 3] vector of coordinates of the best change point
%         errors = [(n-1) x 10] matrix with columns:
%     [ x_cp1 | y_cp1 | z_cp1 | x_cp2 | y_cp2 | z_cp2 | 
%           ...totalSSE_using_that_cp | SSE_first_segment | SSE_second_segment | SSE_third_segment] 

function [cp1, cp2, xyz1, xyz2, xyz3, errors, minError] = twoChangePoints(xyz, hubNode_coordinates, plot_on)

numPoints = size(xyz, 1);

if numPoints < 7
    error(['Cannot fit 2 changepoints to ' num2str(numPoints) ' points. Need at least 7.'])
end

% To make this easier, make sure your hubNode is the first listed node in
% your list.
if xyz(end,:)==hubNode_coordinates
    xyz = flipud(xyz);
elseif xyz(1,:)~=hubNode_coordinates
    error('The hub node chosen does not start or end this vessel.')
end

% First, try placing at most 1 changepoint.
[~, ~, ~, errors, ~] = changePoints(xyz, hubNode_coordinates, 0);
N = size(errors, 1);
errors = [errors(:, 1:3) NaN(N, 3) errors(:, 4:6) NaN(N, 1)];

% Next, try placing two changepoints 
% Place CP1 at the second edge point from the hubNode.
% Find the optimal location of CP2 using the remaining points.
% Then, move CP1 to the next point, then the next, finding the optimal
% location for CP2 each time and adding the change point
% coordinates and SSE's to the 'errors' matrix as you go.
for i=3:size(xyz, 1)-4
    cp1 = xyz(i, :);
    xyz1 = xyz(1:i, :);
    N1 = size(xyz1, 1);
    [~,SSE1,~]=linReg(xyz1, xyz1(1, :), 0);
    xyz2 = xyz(i:end, :);
    [cp2, ~, ~, errors2, minError] = changePoints(xyz2, cp1, 0);
    [~, minErrorRow] = ismember(cp2, errors2(:, 1:3), 'rows');
    if ~isnan(cp2)
        relSSE2 = errors2(minErrorRow, 5);
        relSSE3 = errors2(minErrorRow, 6);
%        errors = [errors; cp1 cp2 round(((SSE1/N1) + minError), 16) round(SSE1/N1, 16) round(relSSE2, 16) round(relSSE3, 16)]
        errors = [errors; cp1 cp2 round(((SSE1/N1) + minError), 16) round(SSE1/N1, 16) round(relSSE2, 16) round(relSSE3, 16)];
    end
end

errors = unique(errors, 'stable', 'rows');

[minError, I] = min(errors(:, 7));
cp1 = errors(I, 1:3);
cp2 = errors(I, 4:6);
[~, cp1_index]=ismember(cp1,xyz,'rows');
[~, cp2_index]=ismember(cp2,xyz,'rows');

if cp1_index ~= 0
    xyz1 = xyz(1:cp1_index, :);
    if cp2_index ~= 0
        xyz2 = xyz(cp1_index:cp2_index, :);
        xyz3 = xyz(cp2_index:end, :);
    else
        xyz2 = xyz(cp1_index:end, :);
        xyz3 = [];
    end
else
    xyz1 = xyz;
    xyz2 = [];
    xyz3 = [];
end
    
    

% close all

% Plot the vessel with vectors defined by the optimal change point
if plot_on == 1
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3))
    hold on
    if cp1_index>1 && cp2_index>cp1_index
        xyz1 = xyz(1:cp1_index, :);
        xyz2 = xyz(cp1_index:cp2_index, :);
        xyz3 = xyz(cp2_index:end, :);

        linReg(xyz1, hubNode_coordinates, plot_on);
        linReg(xyz2, cp1, plot_on);
        linReg(xyz3, cp2, plot_on);
    elseif cp1_index>1
        xyz1 = xyz(1:cp1_index, :);
        xyz2 = xyz(cp1_index:end, :);
        xyz3=0;

        linReg(xyz1, hubNode_coordinates, plot_on);
        linReg(xyz2, cp1, plot_on);
    else
        xyz1 = xyz;
        xyz2 = 0;
        xyz3 = 0;
        linReg(xyz1, cp, plot_on);
    end
end
end

