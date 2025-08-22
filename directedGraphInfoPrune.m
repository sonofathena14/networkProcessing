function [orientation,newNetwork, connection, arcsC2, maxNumDaughters, maxNumParents] = directedGraphInfoPrune(arcs, nodes)

[orientation]=prune_orientation(arcs, nodes);
[newNetwork,arcsC2]=networkGeneratorPrune(arcs, orientation);
[connection, maxNumDaughters, maxNumParents]=connectivity(nodes, newNetwork);

end