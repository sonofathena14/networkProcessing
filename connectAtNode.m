function [new_vessel]=connectAtNode(arcs, v1ID, v2ID, nodeID)

V1=arcs{1,v1ID};
V2=arcs{1,v2ID};
[~, pos1]=ismember(nodeID, V1(1,1:2));
[~, pos2]=ismember(nodeID, V2(1,1:2));

if pos1==pos2
    % Flip second vessel.
    A=flipud(V2(2:end, :));
    t=V2(1,1);
    V2(1,1)=V2(1,2);
    V2(1,2)=t;
    V2=[ V2(1,:); A ];
    [~, pos2]=ismember(nodeID, V2(1,1:2));
end

if pos1==1 && pos2==2
    V2(1,2)=V1(1,2);
    V2 = [V2; V1(3:end, :)];
    new_vessel=V2;
elseif pos2==1 && pos1==2
    V1(1,2)=V2(1,2);
    V1 = [V1; V2(3:end, :)];
    new_vessel=V1;
else
    msg='The degree 2 node is in the same position in both vessels.';
    error(msg)
end
end