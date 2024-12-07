% computing and plotting lrr
% Input: newNetwork matrix
% Outputs: Lrr = vector of length to radius ratios for each vessel,
%          R = vector of radii for each vessel,
%          L = vector of lengths for each vessel,

function [Lrr, R, L]=LRR(newNetwork)

Lrr=[];
R=[];
L=[];

numVes=length(newNetwork); % number of vessels
data = newNetwork;

for i=1:numVes
    R=[R;data(i,4)];
    L=[L;data(i,5)];
    ratio=data(i,5)/data(i,4);
    Lrr=[Lrr;ratio];
end 

mouse=input('Which mouse? ', 's');

figure(1);
scatter(R,L);
title('Length to radius ratio vs. vessel radius')
xlabel('vessel radius') 
ylabel('Lrr','FontSize',16)
filename1 = strcat('Lrr_vs_R_',mouse,'.fig');
savefig(filename1)

figure(2);
scatter(L,Lrr);
title('Length to radius ratio vs. vessel length')
xlabel('vessel length') 
ylabel('Lrr','FontSize',16)
filename2 = strcat('Lrr_vs_L_',mouse,'.fig');
savefig(filename2)

figure(5)
loglog(R,Lrr,'o')
title('Log(Lrr) vs. Log(radius)')
xlabel('log(vessel radius)','FontSize',32) 
ylabel('log(vessel length)','FontSize',32)
filename5 = strcat('LnLrr_vs_LnR_',mouse,'.fig');
savefig(filename5)

figure(6)
scatter(R,L,'o');
title('Length vs. radius')
xlabel('vessel radius','FontSize',32) 
ylabel('vessel length','FontSize',32)
filename6 = strcat('L_vs_R_',mouse,'.fig');
savefig(filename6)

end