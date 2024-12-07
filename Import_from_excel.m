function [details] = Import_from_excel(FILE, details)
% Note that CLFILE is the .mat file, and FILE is the .xlsx file

% Read in Excel file and make into cell
fn          = FILE;        
tBC         = [];
opt         = detectImportOptions(fn);
opt.VariableTypes = {'double','double','double','double','double','double','double','double','double','double'};
shts        = sheetnames(fn);
new_details = cell(length(shts),2);

for j = 1:length(shts)
    new_details{j,1} = [shts(j)];
    new_details{j,2} = [tBC;readtable(fn,opt,'Sheet',shts(j))];
end

details{1,9}  = [];
details{1,9}  = 'CP COORDINATES';
details{1,10} = [];
details{1,10} = 'Coefficients';
 
% Overwrite old vessel_details with new components
for s = 1:height(new_details)
     details{s+1,9}       = zeros(nnz(~isnan(new_details{s,2}.cpx)),3);
     details{s+1,10}      = zeros(nnz(~isnan(new_details{s,2}.c)),1);
     details{s+1,9}(:,1)  = new_details{s,2}.cpx(1:nnz(~isnan(new_details{s,2}.cpx))); %presostium
     details{s+1,9}(:,2)  = new_details{s,2}.cpy(1:nnz(~isnan(new_details{s,2}.cpy))); %presostium
     details{s+1,9}(:,3)  = new_details{s,2}.cpz(1:nnz(~isnan(new_details{s,2}.cpz)));
     details{s+1,10}(:,1) = new_details{s,2}.c(1:nnz(~isnan(new_details{s,2}.c))); %preostium
end

% save(data_name,'vessel_details');