Name = 'm2p1_053007';

load(strcat('Networks/Network_Vessels_',Name,'.mat'), 'nodesC2');

x = nodesC2(:,2);
y = nodesC2(:,3);
z = nodesC2(:,4);

k = convhull(x,y,z);

trisurf(k,x,y,z,'FaceColor','cyan')
axis equal