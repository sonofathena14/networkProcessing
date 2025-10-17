close all
Name = 'm1p1_053107';
arcs = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'arcsC3');
nodes = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'nodesC2');
arcs = arcs.arcsC3;
nodes = nodes.nodesC2;
% ==== USER INPUT SECTION ====
startVessels = [
    3954 3900; %left
    4071 685;% superior
    4179 4304; %middle 1
    4018 4379; %middle 2
    4018 4019; %inferior
    4296 1851 %post caval
];

colors = {'r','g','m','m','c','k','g'};  % Define up to 6 distinct colors
fignb = 2;

legendEntries = {
    'Left', 'Superior', 'Middle', 'Inferior', 'Post Caval'
};

NumVessels = length(arcs);
arcMap = zeros(NumVessels, 2);  % Store from/to for each arc

% Build arcMap from arcs
for i = 1:NumVessels
    arcMap(i, :) = arcs{1,i}(1, 1:2);  % First row of each arc contains [from to]
end

% Assign default color
arcColors = repmat({'b'}, NumVessels, 1);  % Default color is blue

% Color arcs downstream of each starting vessel
for s = 1:min(size(startVessels,1), length(colors))
    fromTo = startVessels(s,:);
    downstream = findDownstreamArcs(arcMap, fromTo);  % Find downstream arcs
    for i = downstream
        arcColors{i} = colors{s};
    end
end

% Plot each arc individually with color
figure(fignb); clf; hold on
plotHandles = gobjects(length(legendEntries),1);  % Preallocate

% Map colors to legend index (note both magenta arcs map to 'Middle')
colorToLegendIndex = containers.Map( ...
    {'r', 'g', 'm', 'c', 'k'}, ...
    {1,   2,   3,   4,   5} ...
);

legendHandles = containers.Map();

for i = 1:NumVessels
    arc = arcs{1,i};
    coords = arc(2:end, 1:3);
    color = arcColors{i};
    
    % Plot arc
    h = plot3(coords(:,1), coords(:,2), coords(:,3), ...
              'Color', color, 'LineWidth', 5);
    
    % Store one handle per relevant color for the legend
    if ~strcmp(color, 'b') && ~isKey(legendHandles, color)
        if isKey(colorToLegendIndex, color)
            idx = colorToLegendIndex(color);
            legendHandles(color) = h;
            plotHandles(idx) = h;
        end
    end
end

% Remove unused handles
validIdx = plotHandles ~= gobjects(1);
plotHandles = plotHandles(validIdx);
legendLabels = legendEntries(validIdx);

legend(plotHandles, legendLabels, 'Location', 'best');

% Plot nodes and labels
% label = num2str(nodes(:,1));
% xN = nodes(:,2); yN = nodes(:,3); zN = nodes(:,4);
% t = text(xN, yN, zN, label, 'FontSize', 16);
% 
% for i = 1:length(nodes)
%     if nodes(i,5)==4
%         t(i).Color = 'red';
%         t(i).FontWeight = 'bold';
%     else
%         t(i).Color = 'black';
%     end
% end

% Extract unique node IDs from the input vessels
uniqueNodes = unique(startVessels(:));  % all 'from' and 'to' nodes flattened and unique

% Find indices of those nodes in the nodes array (assuming nodes(:,1) is node IDs)
[~, loc] = ismember(uniqueNodes, nodes(:,1));

% Coordinates and labels for only these nodes
xN = nodes(loc, 2);
yN = nodes(loc, 3);
zN = nodes(loc, 4);
label = num2str(nodes(loc, 1));

% Plot labels for these nodes only
t = text(xN, yN, zN, label, 'FontSize', 30);

for i = 1:length(t)
    nodeIndex = loc(i);
    if nodes(nodeIndex, 5) == 4
        t(i).Color = 'red';
        t(i).FontWeight = 'bold';
    else
        t(i).Color = 'black';
    end
end

view(0,0)
set(gca, 'fontsize', 30);
hold off

function downstream = findDownstreamArcs(arcMap, startArc)
    downstream = []; 
    visited = false(size(arcMap,1),1);

    startIdx = find(arcMap(:,1) == startArc(1) & arcMap(:,2) == startArc(2), 1);
    if isempty(startIdx)
        warning('Starting arc [%d %d] not found.', startArc(1), startArc(2));
        return;
    end

    queue = startIdx;
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];
        if visited(current)
            continue;
        end
        visited(current) = true;
        downstream(end+1) = current;

        toNode = arcMap(current, 2);
        nextArcs = find(arcMap(:,1) == toNode);
        queue = [queue; nextArcs(:)];
    end
end