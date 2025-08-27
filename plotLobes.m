close all
Name = 'm1p2_060107';
arcs = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'arcsC3');
nodes = load(strcat('Networks/Network_Vessels_',Name,'.mat'),'nodesC2');
arcs = arcs.arcsC3;
nodes = nodes.nodesC2;
% ==== USER INPUT SECTION ====
startVessels = [
    2089 2169; %left
    2162 2163;% superior
    2164 2098; %middle 1
    2153 927; %middle 2
    2153 2151; %inferior
    2164 2190 %post caval
];

colors = {'r','g','m','m','c','k','g'};  % Define up to 6 distinct colors
fignb = 2;

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
for i = 1:NumVessels
    arc = arcs{1,i};
    coords = arc(2:end, 1:3);
    plot3(coords(:,1), coords(:,2), coords(:,3), ...
          'Color', arcColors{i}, 'LineWidth', 1.5);
end

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
t = text(xN, yN, zN, label, 'FontSize', 16);

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
set(gca, 'fontsize', 16);
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