function [DataOut] = get_radii_aorta(vessel_details,Data,Scale)
%function get_radii_aorta; (details, Scale, Data)
%close all;

%load Output/Segmentation2_control_14.mat
DataOut = Data;
details = vessel_details;

disp('in get_radii_aorta')
NVes = length(details(:,1))-1;
for i = 2:NVes+1
  VesNum(i-1) = str2num(details{i,1});
end
lengths = cell2mat(details(2:end,3));

VID(1)  = str2double(details{2,1});
VIDN(1) = find(VesNum==VID(1))+1;

[~, maxIndex] = max(lengths);
VID(3)  = str2double(details{maxIndex + 1, 1});
VIDN(3) = find(VesNum==VID(3))+1;

VID(2)  = details{VIDN(3), 7};
VIDN(2) = find(VesNum==VID(2))+1;

Scale = Scale/10; 
for i = 1:3
   RVes{1,i} = details{VIDN(i),2}(2:end,4)*Scale;
   LVes{1,i} = details{VIDN(i),2}(2:end,5)*Scale;
end

radii_tot  = [RVes{1,1}(:)' RVes{1,2}(2:end)' RVes{1,3}(2:end)']';
xStart     = [0 LVes{1,1}(end) LVes{1,1}(end)+LVes{1,2}(end)];

length_tot = [LVes{1,1}(:)' xStart(2)+LVes{1,2}(2:end)' xStart(3)+LVes{1,3}(2:end)']';

disp('Vessel IDs along the aorta:')
disp([VID(1) VID(2) VID(3)]);
figure(100);
plot(length_tot,radii_tot,'b*','linewidth',3); hold on;
plot(LVes{1,1}(end),RVes{1,1}(end),'mo',LVes{1,1}(end)+LVes{1,2}(end),RVes{1,2}(end),'mo',...
     LVes{1,1}(end)+LVes{1,2}(end)+LVes{1,3}(end),RVes{1,3}(end),'mo','linewidth',3,'markersize',7);
set(gca,'fontsize',16);
xlabel('length (cm)');
ylabel('radius (cm)');
title('Aorta');

TaperIDs = Data.TaperIDs;
if ~isempty(TaperIDs)
    TapID = 2;
else
    TapID = 1;
    TaperIDs = [];
end

vesname = {"ascending aorta","aortic arch","descending aorta"};
for i = 1:3
    prompt = strcat("Should the ",vesname{i}," be split: y(1)/n(0) default?");
    xspl = input(prompt);
    if isempty(xspl)
        xspl = 0;
    end
    if xspl == 1
        prompt = "How many segments (interger, default n=2), click at location(s)?";
        n = input(prompt);
        if isempty(n)
             n = 2;
        end
        [xin,~] = ginput(n-1);
        ID(1) = 1; 
        xStartSeg(1) = 0;
        for j = 1:n
            if j < n
              ID(j+1)= find(LVes{1,i}(:)+xStart(i)<xin(j),1,'last');
            else
              ID(j+1)= length(LVes{1,i}(:));
            end
            Lseg{i,j} = LVes{1,i}(ID(j):ID(j+1))-LVes{1,i}(ID(j));
            Rseg{i,j} = RVes{1,i}(ID(j):ID(j+1));
            
            xStartSeg(j+1) = LVes{1,i}(ID(j+1));
            figure(100); hold on;
            plot(Lseg{i,j}(:)+xStartSeg(j)+xStart(i),Rseg{i,j}(:),'*','linewidth',3); 
            set(gca,'fontsize',16);
            xlabel('length (cm)');
            ylabel('radius (cm)');
            title('Aorta');
            
            prompt = strcat("Should segment ",num2str(j)," of the ",vesname{i}," taper: y(1)/n(0) default?");
            x = input(prompt);
            if x == 1
                xTa(i,j) = 1;
            else
                xTa(i,j) = 0;
            end
            if xTa(i,j) == 0
                Rad = mean(Rseg{i,j}(:));
                figure(100); hold on;
                plot(Lseg{i,j}(:)+xStartSeg(j)+xStart(i),Rad*ones(size(Lseg{i,j}(:))),'r--','linewidth',3);
                set(gca,'fontsize',16);
                xlabel('length (cm)');
                ylabel('radius (cm)');
                title('Aorta');
                pause(1);
                
                prompt = "should the current radius (see Figure) be changed: y(1)/n(0) default?";
                xRad = input(prompt);
                if isempty(xRad)
                    xRad = 0;
                end
                if xRad == 1
                    disp('Click on a radius value for the highlighted segment');
                    [~,Rad] = ginput(1);
                    
                    figure(100); hold on;
                    plot(Lseg{i,j}(:)+xStartSeg(j)+xStart(i),Rad*ones(size(Lseg{i,j})),'r','linewidth',3);
                    set(gca,'fontsize',16);
                    xlabel('length (cm)');
                    ylabel('radius (cm)');
                    title('Aorta');
                end
                Rcomp{i,j}(1) = 0;
                Rcomp{i,j}(2) = 0;
                Rcomp{i,j}(3) = Rad;

                DataOut.rin(ID) = Rad;
                DataOut.rout(ID)= Rad;
                DataOut.error(ID) = std(Rseg{i,j}(:));
                disp([ID DataOut.rin(ID) DataOut.rout(ID)]);
           
            else
                if Rseg{i,j}(1) > Rseg{i,j}(end) 
                    [k1,k2,k3] = find_tapering(Lseg{i,j}(:),Rseg{i,j}(:),xStartSeg(j)+xStart(i));
                    Rcomp{i,j}(1) = k1;
                    Rcomp{i,j}(2) = k2;
                    Rcomp{i,j}(3) = k3;
                else
                    [k1,k2,k3] = find_expanding(Lseg{i,j},Rseg{i,j},xStartSeg(j)+xStart(i));
                    Rcomp{i,j}(1) = k1;
                    Rcomp{i,j}(2) = -k2;
                    Rcomp{i,j}(3) = k3;
                end
            end    
        end % for each segment
        
        if n == 2
            ID   = find(Data.mapIDs(:,2)==VID(i),1,"first");
            dIDs = Data.connectivity(ID,2:4);
            DataOut.connectivity(ID,2)       = NVes;
            DataOut.connectivity(ID,3:4)     = 0;           
            DataOut.connectivity(NVes+1,1)   = NVes;
            DataOut.connectivity(NVes+1,2:4) = dIDs;
            
            tID = find(Data.TermVes==ID-1,1,"first");
            if ~isempty(tID)
                DataOut.TermVes = [DataOut.TermVes(1:tID-1)' DataOut.TermVes(tID+1:end)' NVes]';
            end
            
            DataOut.mapIDs(NVes+1,1) = NVes;
            DataOut.mapIDs(NVes+1,2) = NVes+1;
            DataOut.lengths(ID)      = Lseg{i,1}(end);
            DataOut.lengths(NVes+1)  = Lseg{i,2}(end);
            
            % Radii segment 1
            k1 = Rcomp{i,1}(1);
            k2 = Rcomp{i,1}(2);
            k3 = Rcomp{i,1}(3);
            Rad1 = k1*exp(-k2*Lseg{i,1}(:))+k3;
            DataOut.rin(ID)   = Rad1(1);
            DataOut.rout(ID)  = Rad1(end);
            disp([ID DataOut.rin(ID) DataOut.rout(ID)]);
            DataOut.error(ID) = std(Rseg{i,1});
            if xTa(i,1) ~=0
               TaperIDs(TapID,1) = ID-1;
               TaperIDs(TapID,2) = k1;
               TaperIDs(TapID,3) = k2;
               TaperIDs(TapID,4) = k3;
               TapID = TapID+1;
            end

            % Radii segment 2
            k1 = Rcomp{i,2}(1);
            k2 = Rcomp{i,2}(2);
            k3 = Rcomp{i,2}(3);
            Rad2 = k1*exp(-k2*Lseg{i,2}(:))+k3;
            DataOut.rin(NVes+1)   = Rad2(1);
            DataOut.rout(NVes+1)  = Rad2(end);
            DataOut.error(NVes+1) = std(Rseg{i,2});
            disp([ID DataOut.rin(ID) DataOut.rout(ID)]);
            if xTa(i,2) ~=0
                 TaperIDs(TapID,1) = NVes;
                 TaperIDs(TapID,2) = k1;
                 TaperIDs(TapID,3) = k2;
                 TaperIDs(TapID,4) = k3;
                 TapID = TapID+1;
            end
            NVes = NVes + 1;
        end % Update output  n = 2%

        if n == 3
            ID  = find(Data.mapIDs(:,2)==VID(i),1,"first");
            dIDs = Data.connectivity(ID,2:4);
            DataOut.connectivity(ID,2)       = NVes;
            DataOut.connectivity(ID,3:4)     = 0;   
            DataOut.connectivity(NVes+1,1)   = NVes;
            DataOut.connectivity(NVes+1,2)   = NVes+1;
            DataOut.connectivity(NVes+1,3:4) = 0;
            DataOut.connectivity(NVes+2,1)   = NVes+1;
            DataOut.connectivity(NVes+2,2:4) = dIDs;
            
            tID = find(Data.TermVes==ID-1,1,"first");
            if ~isempty(tID)
                DataOut.TermVes = [DataOut.TermVes(1:tID-1)' DataOut.TermVes(tID+1:end)' NVes+1]';
            end
            
            DataOut.mapIDs(NVes+1,1) = NVes;
            DataOut.mapIDs(NVes+2,1) = NVes+1;
            DataOut.mapIDs(NVes+1,2) = NVes+1;
            DataOut.mapIDs(NVes+2,2) = NVes+2;
            DataOut.lengths(ID)      = Lseg{i,1}(end);
            DataOut.lengths(NVes+1)  = Lseg{i,2}(end);
            DataOut.lengths(NVes+2)  = Lseg{i,3}(end);
  
            % Radii segment 1
            k1 = Rcomp{i,1}(1);
            k2 = Rcomp{i,1}(2);
            k3 = Rcomp{i,1}(3);
            Rad1 = k1*exp(-k2*Lseg{i,1}(:))+k3;
            DataOut.rin(ID)   = Rad1(1);
            DataOut.rout(ID)  = Rad1(end);
            DataOut.error(ID) = std(Rseg{i,1});
            disp([ID DataOut.rin(ID) DataOut.rout(ID)]);
            if xTa(i,1) ~=0
               TaperIDs(TapID,1) = ID-1;
               TaperIDs(TapID,2) = k1;
               TaperIDs(TapID,3) = k2;
               TaperIDs(TapID,4) = k3;
               TapID = TapID+1;
            end
           
            % Radii segment 2
            k1 = Rcomp{i,2}(1);
            k2 = Rcomp{i,2}(2);
            k3 = Rcomp{i,2}(3);
            Rad2 = k1*exp(-k2*Lseg{i,2}(:))+k3;
            DataOut.rin(NVes+1)   = Rad2(1);
            DataOut.rout(NVes+1)  = Rad2(end);
            DataOut.error(NVes+1) = std(Rseg{i,2});
            disp([NVes+1 DataOut.rin(NVes+1) DataOut.rout(NVes+1)])
            
            if xTa(i,2) ~=0
                 TaperIDs(TapID,1) = NVes;
                 TaperIDs(TapID,2) = k1;
                 TaperIDs(TapID,3) = k2;
                 TaperIDs(TapID,4) = k3;
                 TapID = TapID+1;
            end

            % Radii segment 3
            k1 = Rcomp{i,3}(1);
            k2 = Rcomp{i,3}(2);
            k3 = Rcomp{i,3}(3);
            Rad3 = k1*exp(-k2*Lseg{i,3}(:))+k3;
            DataOut.rin(NVes+2)   = Rad3(1);
            DataOut.rout(NVes+2)  = Rad3(end);
            DataOut.error(NVes+2) = std(Rseg{i,3});
            disp([NVes+2 DataOut.rin(NVes+2) DataOut.rout(NVes+2)]);
            if xTa(i,3) ~=0
                 TaperIDs(TapID,1) = NVes+1;
                 TaperIDs(TapID,2) = k1;
                 TaperIDs(TapID,3) = k2;
                 TaperIDs(TapID,4) = k3;
                 TapID = TapID+1;
            end
        end
            
    else % no splitting
        prompt = strcat("Should ",vesname{i}," taper: y(1)/n(0) default?");
        x = input(prompt);
        if x == 1
            xTa(i,1) = 1;
        else
            xTa(i,1) = 0;
        end
        if xTa(i,1) == 0
            Rad = mean(RVes{1,i}(:));
            figure(100); hold on;
            plot(LVes{1,i}(:)+xStart(i),Rad*ones(size(LVes{1,i}(:))),'r--','linewidth',3);
            set(gca,'fontsize',16);
            xlabel('length (cm)');
            ylabel('radius (cm)');
            title('Aorta');
            pause(1);
            
            prompt = "should the current radius (see Figure) be changed: y(1)/n(0) default?";
            xRad = input(prompt);
            if isempty(xTa(i,1))
                xRad = 0;
            end
            if xRad == 1
            disp('Click on a radius value for the highlighted segment');
            [~,Rad] = ginput(1);
                
            figure(100); hold on;
            plot(LVes{1,i}(:)+xStart(i),Rad*ones(size(LVes{1,i})),'r','linewidth',3);
            set(gca,'fontsize',16);
            xlabel('length (cm)');
            ylabel('radius (cm)');
            title('Aorta');
            
            Rcomp{i,1}(1) = 0;
            Rcomp{i,1}(2) = 0;
            Rcomp{i,1}(3) = Rad;
            
            ID  = find(Data.mapIDs(:,2)==VID(i),1,"first");
            DataOut.rin(ID)   = Rad;
            DataOut.rout(ID)  = Rad;
            disp([ID DataOut.rin(ID) DataOut.rout(ID)]);
            DataOut.error(ID) = std(RVes{1,i}(:));
            
            end % if xRad == 1
        else 
            if RVes{1,i}(1) > RVes{1,i}(end) 
                [k1,k2,k3] = find_tapering(LVes{1,i}(:),RVes{1,i}(:),xStart(i));
                Rcomp{i,1}(1) = k1;
                Rcomp{i,1}(2) = k2;
                Rcomp{i,1}(3) = k3;
            else
                [k1,k2,k3] = find_expanding(LVes{1,i}(:),RVes{1,i}(:),xStart(i));
                Rcomp{i,1}(1) = k1;
                Rcomp{i,1}(2) = -k2;
                Rcomp{i,1}(3) = k3;
            end % end if

            ID  = find(Data.mapIDs(:,2)==VID(i),1,"first");
            k1 = Rcomp{i,1}(1);
            k2 = Rcomp{i,1}(2);
            k3 = Rcomp{i,1}(3);
            Rad = k1*exp(-k2*LVes{1,i}(:))+k3;
            DataOut.rin(ID)   = Rad(1);
            DataOut.rout(ID)  = Rad(end);
            DataOut.error(ID) = std(RVes{1,i});
            disp([ID DataOut.rin(ID) DataOut.rout(ID)]);
            TaperIDs(TapID,1)   = ID-1;
            TaperIDs(TapID,2)   = k1;
            TaperIDs(TapID,3)   = k2;
            TaperIDs(TapID,4)   = k3;
            TapID = TapID+1;
        end % if tapering vessel j
    end % vessel i
end % for each vessel

if ~isempty(TaperIDs)
  DataOut.TaperIDs = TaperIDs;
end
