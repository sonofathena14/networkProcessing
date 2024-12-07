function [pressure, flow] = ExtractData(nv,ntp,file_id)

% NOTE: Results are stored such that the first 1:N entries are the PROXIMAL
% large artery (1-15) and large vein (16-27) predictions, N+1:2N are the
% MIDPOINT predictions in the same vessels, and 2N+1:3N are the DISTAL
% predictions.

data = load(sprintf('output_%d.2d', file_id)); 

[~,~,p,q,~,~] = gnuplot(data);

pressure_all = (p(:,nv+1:2*nv)); % middle prediction
flow_all     = (q(:,nv+1:2*nv)); % middle prediction
flow_steps   = round(linspace(1,size(q,1),ntp));
flow_all     = flow_all(flow_steps,:);

% real data consists of flow (20 equispaced points) from vessels 2-7, 
% and pressure (min and max points only) from vessel 9, 
% so these are also the predictions we save
pressure  = [max(pressure_all(:,9)), min(pressure_all(:,9))]; %???SHOULD THIS BE (:,9) OR (9,:)?
flow = flow_all(:,[3,5:7]); % skip inlet vessel
maxPredFlow = [max(flow(:,1)), max(flow(:,2)), max(flow(:,3)), max(flow(:,4))];
minPredFlow = [min(flow(:,1)), min(flow(:,2)), min(flow(:,3)), min(flow(:,4))];

end