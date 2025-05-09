function volumes = edgeVolume(vessel_details)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

volumes = zeros(length(vessel_details)-1,1);
for i=2:(length(vessel_details))
    edgeVol = 0;
    %disp(i)
    %disp(length(arcsC3{1,i})-1)
    for j=2:(size(vessel_details{i,2},1)-1)
        %disp(j)
        bot = vessel_details{i,2}(j,4);
        top = vessel_details{i,2}(j+1,4);
        len = vessel_details{i,2}(j+1,5)- vessel_details{i,2}(j,5);
        conVol = pi*len/3*(top^2+bot^2+top*bot);
        edgeVol = edgeVol+conVol;
    end
    volumes(i-1,1) = edgeVol;
end

end