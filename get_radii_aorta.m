function DataOut = get_radii_aorta(details, Scale, DataFC)
close all;

global data;

DataOut = DataFC;
NVes    = length(details(:,1));
for i = 2:NVes
  VesNum(i-1) = str2num(details{i,1});
end

AsID = str2double(details{2,1});
AsIDN = find(VesNum==AsID)+1;

lengths = cell2mat(details(2:end,3));  
[~, maxIndex] = max(lengths);
DesID = str2double(details{maxIndex + 1, 1});
DesIDN = find(VesNum==DesID)+1;

ArID  = details{DesIDN, 7};
ArIDN  = find(VesNum==ArID)+1;

Scale = Scale/10; 
radii_tot  = [details{AsIDN,2}(2:end,4)' details{ArIDN,2}(2:end,4)' details{DesIDN,2}(2:end,4)']';
length_tot = [details{AsIDN,2}(2:end,5)' details{ArIDN,2}(2:end,5)'+details{AsIDN,2}(end,5) details{DesIDN,2}(2:end,5)'+details{AsIDN,2}(end,5)+details{ArIDN,2}(end,5)]';
radii_tot  = radii_tot*Scale;
length_tot = length_tot*Scale;

RVes1   = details{AsIDN,2}(2:end,4)*Scale;
LVes1   = details{AsIDN,2}(2:end,5)*Scale;
NVes1   = length(LVes1);

RVes2   = details{ArIDN,2}(2:end,4)*Scale;
LVes2   = details{ArIDN,2}(2:end,5)*Scale;
NVes2   = length(LVes2);

RVes3   = details{DesIDN,2}(2:end,4)*Scale;
LVes3   = details{DesIDN,2}(2:end,5)*Scale;
NVes3   = length(LVes3);

figure(100);
plot(length_tot,radii_tot,'b*','linewidth',3,'markersize',7); hold on;
plot(LVes1(end),details{AsIDN,2}(end,4)*Scale,'mo',LVes1(end)+LVes2(end),details{ArIDN,2}(end,4)*Scale,'mo','linewidth',4,'markersize',10 );
set(gca,'fontsize',16);
xlabel('length (cm)');
ylabel('radius (cm)');
title('Aorta');

opts = optimset('MaxIter',10000,'MaxFunEvals',10000);
k2I  = 0.1;

prompt = "Should ascending aorta taper (y(1)/n(0))?";
xAs = input(prompt);
iTapVes = 0;

if xAs == 1
  disp('Click on Rin value');
  [x,~] = ginput(1);

  ID = max(find(LVes1<x))+1;
  data.x = LVes1;
  data.r = [RVes1(ID)*ones(size(RVes1(1:ID-1)))' RVes1(ID:end)']';
  
  opts = optimset('MaxIter',10000,'MaxFunEvals',10000);
  [k2As, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
  cstAs  =  exp(-k2As*data.x(end));
  k3As   = (data.r(end)-data.r(1)*cstAs)/(1-cstAs);
  k1As   = data.r(1)-k3As;
 
  RestVes1 = k1As*exp(-k2As*LVes1)+k3As;
  
  N1 = length(LVes1);
  sumx  = sum((data.r-RestVes1).^2);
  error = sqrt(1/(N1-1)*sumx);

  figure(100); hold on;
  plot(LVes1,RestVes1,'m','linewidth',4,'markersize',10);
  set(gca,'fontsize',16);
  xlabel('length (cm)');
  ylabel('radius (cm)');
  title('Aorta');
  
  rID = find(DataFC.mapIDs(:,2)==AsID);
  nID = DataFC.mapIDs(rID,1);
  DataOut.rin(nID+1) = RestVes1(1);
  DataOut.rout(nID+1)= RestVes1(end);
  DataOut.error(nID+1) = error;
 
  TaperIDs(iTapVes+1,1) = nID;
  TaperIDs(iTapVes+1,2) = k1As;
  TaperIDs(iTapVes+1,3) = k2As;
  TaperIDs(iTapVes+1,4) = k3As;
  iTapVes = iTapVes+1;
end

prompt = "Should the aortic arch taper (y(1)/n(0))?";
xAr = input(prompt);
if xAr == 1
  data.x = LVes2;
  data.r = RVes2;

  [k2Ar, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
  cstAr  =  exp(-k2Ar*data.x(end));
  k3Ar   = (data.r(end)-data.r(1)*cstAr)/(1-cstAr);
  k1Ar   = data.r(1)-k3Ar;

  RestVes2 = k1Ar*exp(-k2Ar*LVes2)+k3Ar;

  N2 = length(LVes2);
  sumx  = sum((data.r-RestVes2).^2);
  error = sqrt(1/(N2-1)*sumx);

  figure(100); 
  plot(LVes2+LVes1(end),RestVes2,'m','linewidth',4,'markersize',10);
  set(gca,'fontsize',16);
  xlabel('length (cm)');
  ylabel('radius (cm)');
  title('Aorta');

  rID = find(DataFC.mapIDs(:,2)==ArID);
  nID = DataFC.mapIDs(rID,1);
  DataOut.rin(nID+1)   = RestVes2(1);
  DataOut.rout(nID+1)  = RestVes2(end);
  DataOut.error(nID+1) = error;
  
  TaperIDs(iTapVes+1,1) = nID;
  TaperIDs(iTapVes+1,2) = k1Ar;
  TaperIDs(iTapVes+1,3) = k2Ar;
  TaperIDs(iTapVes+1,4) = k3Ar;
  iTapVes = iTapVes+1;
end

prompt = "Should the descending aorta be split (y(1)/n(0))?";
xspl = input(prompt);
if xspl == 1
    prompt = "How many segments is needed (either 2 or 3)";
    n = input(prompt);
    [x,~] = ginput(n-1);
  
    if n == 2
      ID = max(find(LVes3<x(1)-LVes1(end)-LVes2(end)));
      Lseg1 = LVes3(1:ID);
      Lseg2 = LVes3(ID:end)-LVes3(ID);
      Rseg1 = RVes3(1:ID);
      Rseg2 = RVes3(ID:end);

      data.x = Lseg1;
      data.r = Rseg1;
      [k2s1, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
      csts1  =  exp(-k2s1*data.x(end));
      k3s1   = (data.r(end)-data.r(1)*csts1)/(1-csts1);
      k1s1   = data.r(1)-k3s1;
      RestSeg1 = k1s1*exp(-k2s1*Lseg1)+k3s1;

      Ns1 = length(Rseg1);
      sumx1  = sum((data.r-RestSeg1).^2);
      errors1 = sqrt(1/(Ns1-1)*sumx1);

      disp('Click on Rin value');
      figure(100); 
      [x,~] = ginput(1);
      ID = max(find(Lseg2<x-LVes1(end)-LVes2(end)-Lseg1(end)))+1;
      
      data.x = Lseg2;
      data.r = [Rseg2(1:ID)' Rseg2(ID)*ones(size(Rseg2(ID+1:end)))']';
      [k2s2, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
      csts2  =  exp(-k2s2*data.x(end));
      k3s2   = (data.r(end)-data.r(1)*csts2)/(1-csts2);
      k1s2   = data.r(1)-k3s2;
      RestSeg2 = k1s2*exp(-k2s2*Lseg2)+k3s2;
      
      Ns2 = length(Rseg2);
      sumx2  = sum((data.r-RestSeg2).^2);
      errors2 = sqrt(1/(Ns2-1)*sumx2);

      figure(100); 
      plot(Lseg1+LVes1(end)+LVes2(end),RestSeg1,'m','linewidth',4,'markersize',10); 
      plot(Lseg2+LVes1(end)+LVes2(end)+Lseg1(end),RestSeg2,'c','linewidth',4,'markersize',10); 
      set(gca,'fontsize',16);
      xlabel('length (cm)');
      ylabel('radius (cm)');
      title('Aorta');

      % Updating connectivity
      rID = find(DataFC.mapIDs(:,2)==DesID);
      DataOut.connectivity(rID,2)    = NVes-1;
      DataOut.connectivity(NVes,1)   = NVes-1;
      DataOut.connectivity(NVes,2:4) = 0;

      % Updating output segment 1
      DataOut.lengths(rID) = Lseg1(end);
      DataOut.rin(rID)  = RestSeg1(1);
      DataOut.rout(rID) = RestSeg1(end);
      DataOut.error(rID)= errors1;
      
      TaperIDs(iTapVes+1,1) = DesID;
      TaperIDs(iTapVes+1,2) = k1s1;
      TaperIDs(iTapVes+1,3) = k2s1;
      TaperIDs(iTapVes+1,4) = k3s1;
      iTapVes = iTapVes+1;

      % Updating output segment 2
      DataOut.lengths(NVes)= Lseg2(end);
      DataOut.rin(NVes)    = RestSeg2(1);
      DataOut.rout(NVes)   = RestSeg2(end);
      DataOut.error(NVes)  = errors2;
      
      TaperIDs(iTapVes+1,1) = NVes;
      TaperIDs(iTapVes+1,2) = k1s2;
      TaperIDs(iTapVes+1,3) = k2s2;
      TaperIDs(iTapVes+1,4) = k3s2;
      iTapVes = iTapVes+1;

      % Updating terminal vessels
      tID = find(DataOut.TermVes==DataOut.mapIDs(rID,1));
      DataOut.TermVes = [DataOut.TermVes(1:tID-1)' DataOut.TermVes(tID+1:end)' NVes]';

      DataOut.mapIDs(NVes,1)   = NVes-1;
      DataOut.mapIDs(NVes,2)   = NVes;
    elseif n ==3
      ID1 = max(find(LVes3<x(1)-LVes1(end)-LVes2(end)));
      ID2 = max(find(LVes3<x(2)-LVes1(end)-LVes2(end)));
      Lseg1 = LVes3(1:ID1);
      Lseg2 = LVes3(ID1:ID2)-LVes3(ID1);
      Lseg3 = LVes3(ID2:end)-LVes3(ID2);
      Rseg1 = RVes3(1:ID1);
      Rseg2 = RVes3(ID1:ID2);
      Rseg3 = RVes3(ID2:end);
      
      data.x = Lseg1;
      data.r = Rseg1;
      [k2s1, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
      csts1  =  exp(-k2s1*data.x(end));
      k3s1   = (data.r(end)-data.r(1)*csts1)/(1-csts1);
      k1s1   = data.r(1)-k3s1;
      RestSeg1 = k1s1*exp(-k2s1*Lseg1)+k3s1;

      Ns1 = length(Rseg1);
      sumx1  = sum((data.r-RestSeg1).^2);
      errors1 = sqrt(1/(Ns1-1)*sumx1);

      data.x = Lseg2;
      data.r = Rseg2;
      [k2s2, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
      csts2  =  exp(-k2s2*data.x(end));
      k3s2   = (data.r(end)-data.r(1)*csts2)/(1-csts2);
      k1s2   = data.r(1)-k3s2;
      RestSeg2 = k1s2*exp(-k2s2*Lseg2)+k3s2;

      Ns2 = length(Rseg2);
      sumx2  = sum((data.r-RestSeg2).^2);
      errors2 = sqrt(1/(Ns2-1)*sumx2);

      disp('Click on Rin value')
      figure(100);
      [x,~] = ginput(1);
      ID = max(find(Lseg3<x-LVes1(end)-LVes2(end)-Lseg1(end)-Lseg2(end)))+1;
      
      data.x = Lseg3;
      data.r = [Rseg3(1:ID)' Rseg3(ID)*ones(size(Rseg3(ID+1:end)))']';
      [k2s3, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
      csts3  =  exp(-k2s3*data.x(end));
      k3s3   = (data.r(end)-data.r(1)*csts3)/(1-csts3);
      k1s3   = data.r(1)-k3s3;
      RestSeg3 = k1s3*exp(-k2s3*Lseg3)+k3s3;

      Ns3 = length(Rseg3);
      sumx3  = sum((data.r-RestSeg3).^2);
      errors3 = sqrt(1/(Ns3-1)*sumx3);

      figure(100); 
      plot(Lseg1+LVes1(end)+LVes2(end),RestSeg1,'m','linewidth',4,'markersize',10); 
      plot(Lseg2+LVes1(end)+LVes2(end)+Lseg1(end),RestSeg2,'c','linewidth',4,'markersize',10);
      plot(Lseg3+LVes1(end)+LVes2(end)+Lseg1(end)+Lseg2(end),RestSeg3,'g','linewidth',4,'markersize',10);  
      set(gca,'fontsize',16);
      xlabel('length (cm)');
      ylabel('radius (cm)');
      title('Aorta');

      rID = find(DataFC.mapIDs(:,2)==DesID);
      % Updating connectivity matrix
      DataOut.connectivity(rID,2)     = NVes-1;
      DataOut.connectivity(NVes,1)    = NVes-1;
      DataOut.connectivity(NVes+1,1)  = NVes;
      DataOut.connectivity(NVes,2)    = NVes;
      DataOut.connectivity(NVes,3:4)  = 0;
      DataOut.connectivity(NVes+1,2:4)= 0;

      % Updating outputs segment 1
      DataOut.lengths(rID) = Lseg1(end);
      DataOut.rin(rID)  = RestSeg1(1);
      DataOut.rout(rID) = RestSeg1(end);
      DataOut.error(rID)= errors1;
      
      TaperIDs(iTapVes+1,1) = DesID;
      TaperIDs(iTapVes+1,2) = k1s1;
      TaperIDs(iTapVes+1,3) = k2s1;
      TaperIDs(iTapVes+1,4) = k3s1;
      iTapVes = iTapVes+1;
      
      % Updating outputs segment 2
      DataOut.lengths(NVes)= Lseg2(end);
      DataOut.rin(NVes)    = RestSeg2(1);
      DataOut.rout(NVes)   = RestSeg2(end);
      DataOut.error(NVes)  = errors2;
      
      TaperIDs(iTapVes+1,1) = NVes;
      TaperIDs(iTapVes+1,2) = k1s2;
      TaperIDs(iTapVes+1,3) = k2s2;
      TaperIDs(iTapVes+1,4) = k3s2;
      iTapVes = iTapVes+1;

      % Updating outputs segment 3
      DataOut.lengths(NVes+1)= Lseg3(end);
      DataOut.rin(NVes+1)    = RestSeg3(1);
      DataOut.rout(NVes+1)   = RestSeg3(end);
      DataOut.error(NVes+1)  = errors3;

      TaperIDs(iTapVes+1,1) = NVes+1;
      TaperIDs(iTapVes+1,2) = k1s3;
      TaperIDs(iTapVes+1,3) = k2s3;
      TaperIDs(iTapVes+1,4) = k3s3;
      iTapVes = iTapVes+1;

      tID = find(DataOut.TermVes==DataOut.mapIDs(rID,1));
      DataOut.TermVes = [DataOut.TermVes(1:tID-1)' DataOut.TermVes(tID+1:end)' NVes+1]';

      DataOut.mapIDs(NVes,1)   = NVes-1;
      DataOut.mapIDs(NVes,2)   = NVes;
      DataOut.mapIDs(NVes+1,1) = NVes;
      DataOut.mapIDs(NVes+1,2) = NVes+1;
    end 
else
    prompt = "Should the descending aorta taper (y(1)/n(0))?";
    xdes = input(prompt);
    if xdes == 1
       disp('Click on Rin and Rout value')

       figure(100);
       [x,~] = ginput(2);
       IDin  = max(find(LVes3<x(1)-LVes1(end)-LVes2(end)))+1;
       IDout = max(find(LVes3<x(2)-LVes1(end)-LVes2(end)))+1;
        
       data.x = LVes3;
       data.r = [RVes3(IDin)*ones(size(RVes3(1:IDin)))' RVes3(IDin+1:IDout)' RVes3(IDout)*ones(size(RVes3(IDout+1:end)))']';
   
       [k2Des, ~, ~, ~] = fminsearch(@model_fmin_k13fix,k2I,opts);
       cstDes  =  exp(-k2Des*data.x(end));
       k3Des   = (data.r(end)-data.r(1)*cstDes)/(1-cstDes);
       k1Des   = data.r(1)-k3Des;
       RestVes3 = k1Des*exp(-k2Des*LVes3)+k3Des;

       rID = find(DataFC.mapIDs(:,2)==DesID);
       DataOut.rin(rID)=RestVes3(1);
       DataOut.rin(rID)=RestVes3(end);

       TaperIDs(iTapVes+1,1) = DesID;
       TaperIDs(iTapVes+1,2) = k1Des;
       TaperIDs(iTapVes+1,3) = k2Des;
       TaperIDs(iTapVes+1,4) = k3Des;
       iTapVes = iTapVes+1;
       
       figure(100); 
       plot(LVes3+LVes1(end)+LVes2(end),RestVes3,'m','linewidth',4,'markersize',10);
       set(gca,'fontsize',16);
       xlabel('length (cm)');
       ylabel('radius (cm)');
       title('Aorta');
    end
end

if xAs ~= 0 || xAr~= 0 || xspl ~= 0 || xdes ~= 0 
  TaperOrgIDs = DataFC.TaperIDs;
  if ~isempty(TaperIDs) && ~isempty(TaperOrgIDs)
    sameIDs=find(ismember(TaperOrgIDs(:,1),TaperIDs(:,1)));
    allIDs = [1:size(TaperOrgIDs,1)]';
    IDs    = setdiff(allIDs,sameIDs);
    TaperOrgIDs = TaperOrgIDs([IDs],:);
   
    TaperMerged = [TaperIDs; TaperOrgIDs];
    TaperMerged = sortrows(TaperMerged,1);
    DataOut.TaperIDs = TaperMerged;
    TaperMerged = sortrows(TaperMerged,1);
  elseif isempty(TaperOrgIDs)
      TaperMerged = TaperIDs;
      DataOut.TaperIDs = TaperMerged;
  elseif isempty(TaperIDs)
      TaperMerged = TaperOrgIDs;
      DataOut.TaperIDs = TaperMerged;
  else  
      DataOut.TaperIDs = TaperMerged;
  end

end % function %
