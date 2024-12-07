function data = CreateFluidsCodeInput(vessel_radii,details,maxDaughters)
    
    NumVes = length(vessel_radii(2:end,1));

    lengths  = cell2mat(vessel_radii(2:end, 2)); 
    rin      = cell2mat(vessel_radii(2:end, 3));  
    rout     = cell2mat(vessel_radii(2:end, 4));     
    error    = cell2mat(vessel_radii(2:end, 5));
    TaperIDs = cell2mat(vessel_radii(2:end, 7));
    
   [newConnectivity,mapIDs,TermVes] = CreateConnect(details,maxDaughters);
   
    if size(TaperIDs,1) ~= 0 
      for i = 1:size(TaperIDs,1)
         IDs(i) = find(mapIDs(:,2)==TaperIDs(i,1));
      end
      TaperIDs(:,1) = mapIDs(IDs,1);
    end

    TermVes
    
    data.connectivity = newConnectivity;
    data.mapIDs   = mapIDs;
    data.TermVes  = TermVes;
    data.lengths  = lengths;
    data.rin      = rin;
    data.rout     = rout;
    data.error    = error;
    data.TaperIDs = TaperIDs;
end
