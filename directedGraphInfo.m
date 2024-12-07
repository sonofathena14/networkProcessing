function [orientation, newNetwork, connection, arcsC2, maxNumDaughters, maxNumParents] = directedGraphInfo(arcs, nodes, path)

[orientation]=edge_orientation(arcs, path);
[newNetwork,arcsC2]=networkGenerator(arcs, orientation);
[connection, maxNumDaughters, maxNumParents]=connectivity(nodes, newNetwork);

end