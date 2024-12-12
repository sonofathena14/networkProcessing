% computing and plotting alphas and betas
% Input: newNetwork matrix and nodes matrix and scaling factor for radius
% Outputs: A = vector of alpha values for each bifurcation,
%          B = vector of beta values for each bifurcation,
%          R = vector of parent radii for each bifurcation,
%          node_details = table detailing the number of nodes with certain number of daughters
%          alpha = average alpha for whole network 
%          beta = average beta for whole network
%          lrr  = length to radius ratio

function [alphaVec, betaVec, rad_alphabetaVec, lrrVec, rad_lrr, node_details, alpha, beta, lrr, Connection]=alpha_beta(newNetwork, nodes)
A=[];
B=[];
R=[];
Nodes_With_0_Daughters=0;
Nodes_With_3_Daughters=0;
Nodes_With_2_Daughters=0;
Nodes_With_1_Daughters=0;
More_than_3=0;

[Connection, maxD, ~]=connectivity(nodes, newNetwork);

numNodes=length(Connection);
numVes  = length(newNetwork);
fp=8; % 1st parent column

Connection( all(~Connection,2), : ) = [];

for i=1:numNodes
    nodeID = Connection(i,1);
    parent=Connection(i,fp);
    daughters=zeros(maxD,1);
    for j=1:maxD
        daughters(j,1)=Connection(i,3+j);
    end
    %disp([nodeID parent daughters']);
    daughters = daughters(daughters ~= 0);
    currentNumDaughters=length(daughters);
    if currentNumDaughters==2 && parent~=0
        Parent_Ind=find(newNetwork(:,1)==parent);
        rP = newNetwork(Parent_Ind,4);
        D1_Ind = find(newNetwork(:,1)==daughters(1));
        D2_Ind = find(newNetwork(:,1)==daughters(2));
        rD1=newNetwork(D1_Ind,4);
        rD2=newNetwork(D2_Ind,4);
        if rD1<rD2
            t=rD1;
            rD1=rD2;
            rD2=t;
        end
        alpha=rD1/rP;
        beta=rD2/rP;
        A=[A;alpha];
        B=[B;beta];
        R=[R;rP];
        Nodes_With_2_Daughters=Nodes_With_2_Daughters+1;
    elseif parent==0
        message=['Node ',num2str(i-1),' is root node.'];
        disp(message);
        Root_Node = 1;
    elseif currentNumDaughters==3
        message=['Node ',num2str(i-1),' has 3 daughters.'];
        disp(message);
        Nodes_With_3_Daughters=Nodes_With_3_Daughters+1;
    elseif currentNumDaughters==0
        message=['Node ',num2str(i-1),' is a terminal node.'];
        disp(message);
        Nodes_With_0_Daughters=Nodes_With_0_Daughters+1;
    elseif currentNumDaughters==1
        message=['Node ',num2str(i-1),' has only 1 daughter.'];
        disp(message);
        Nodes_With_1_Daughters=Nodes_With_1_Daughters+1;
    else
        message=['Node ',num2str(i-1),' has more than 3 daughters.'];
        disp(message);  
        More_than_3=More_than_3+1;
    end
end

for i = 1:numVes
    vesID = newNetwork(i,1);
    L(i)     = newNetwork(i,5);
    rad(i)   = newNetwork(i,4);
    length_to_radius(i)= L(i)/rad(i);
end

node_details=table(Root_Node, Nodes_With_0_Daughters, Nodes_With_1_Daughters, Nodes_With_2_Daughters, Nodes_With_3_Daughters,More_than_3);

alphaVec = A;
betaVec  = B;
rad_alphabetaVec = R;
lrrVec   = length_to_radius;
radVec   = rad;

alpha = mean(A);
beta  = mean(B);
lrr   = mean(lrrVec);

figure(11);
ax1 = subplot(2,1,1);
scatter(ax1,R,A)
title('Scaling parameter \alpha vs. parent vessel radius','FontSize',16)
ylabel('\alpha','FontSize',16) 
ax2 = subplot(2,1,2);
scatter(ax2,R,B)
title('Scaling parameter \beta vs. parent vessel radius','FontSize',16)
ylabel('\beta','FontSize',16) 
xlabel('Parent vessel radius','FontSize',16)

figure(12);
scatter(R,A,'mo');
hold on
scatter(R,B,'go');
title('Scaling parameters \alpha and \beta vs. parent vessel radius','FontSize',18)
xlabel('Parent radius (vx)','FontSize',16) 
legend({'\alpha','\beta'},'Location','northeast','FontSize',16)
hold off

figure(13);
scatter(radVec,lrrVec,'bo');
title('Length-to-radius','FontSize',18)
xlabel('Vessel radius (vx)','FontSize',16) 
ylabel('lrr (vx)','FontSize',16) 

end