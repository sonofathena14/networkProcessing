% This function generates a list of vessels in a 3D structured tree. This
% structured tree will start with a root node at the origin.

% INPUT: numGen=the number of generations you want your structured tree to
% go down to. *Biggest tree I've been patient enough to wait for is 10, as
% it starts to take quite a while after that. Also the root vessel
% represents generation 0, and the first 2 daughter vessels represent
% generation 1, etc.

% OUTPUTS: Vessels={1 x (numGen+1)} cell array, with each cell {1,i} 
% containing information about the (i-1)th generation. 
% Each cell {1,i} in Vessels is a {2 x 2^(i-1)} cell array, with each cell 
% containing information about a different vessel in that generation. 
% The first row of cells gives the start and end coordinates of each of 
% the 2^(i-1) vessels of the (i-1)th generation. The second row of cells
% gives the [radius; diameter] of the vessel described in the above row in
% that column.
% theta=the out of plane angle theta used to generate that particular tree,
% taken from a uniform distribution from [0,pi]. 

%generate_3D_tree_MAIN(10)

function [Vessels, theta]=generate_3D_tree(numGen)

tic
%% Initialize variables
d_root=1; % Diameter of the root vessel
l_root=25; % Length of the root vessel
alpha=0.87;
beta=0.69;
theta=random('Uniform',0,pi);

%% Define array to contain vessel start and end points, diameter, length per generation
Vessels=cell(1,numGen+1);
root_start=[0 0 0];
root_end=[0 0 l_root];
Vessels{1,1}=cell(2,1);

%% Populate Vessels array for root
Vessels{1,1}{1,1}=[root_start; root_end];
Vessels{1,1}{2,1}=[d_root;l_root];

%% For each generation:
for i=2:numGen+1
    branches_this_generation=2^(i-1);
    Vessels{1,i}=cell(2,branches_this_generation);
    
    %% For each branch: First populate starting pts, diameters, & lengths
    for j=1:branches_this_generation
        
        if rem(j,2)==1 % odd, alpha
            p_i=(j+1)/2; % Parent index
            s_f=alpha; % Scaling factor
        else % even, beta
            p_i=j/2; % Parent index
            s_f=beta; % Scaling factor
        end
        
        % Set up points cell, record starting point
        Vessels{1,i}{1,j}=zeros(2,3);
        Vessels{1,i}{1,j}(1,:)=Vessels{1,i-1}{1,p_i}(2,:);
        
        % Current branch diameter and length
        Vessels{1,i}{2,j}=s_f*Vessels{1,i-1}{2,p_i};
    end
    
    %% For each branch: Next calculate ending point of each branch
    for j=1:branches_this_generation
        
        % Determine R_psi and current_length
        if rem(j,2)==1 % odd, d1
            p_i=(j+1)/2; % Parent index
            s_i=j+1; % Sibling index
            d0=Vessels{1,i-1}{2,p_i}(1,1);
            d1=Vessels{1,i}{2,j}(1,1);
            d2=Vessels{1,i}{2,s_i}(1,1);
            
            psi=acos(((d0^4)+(d1^4)-(d2^4))/(2*(d0^2)*(d1^2)));
            current_length=Vessels{1,i}{2,j}(2,1);
            
        else % even, d2
            p_i=j/2; % Parent index
            s_i=j-1; % Sibling index
            d0=Vessels{1,i-1}{2,p_i}(1,1);
            d1=Vessels{1,i}{2,s_i}(1,1);
            d2=Vessels{1,i}{2,j}(1,1);
            
            psi=-acos(((d0^4)+(d2^4)-(d1^4))/(2*(d0^2)*(d2^2)));
            current_length=Vessels{1,i}{2,j}(2,1);
        end
        
        x0=Vessels{1,i}{1,j}(1,:)';
        x0_prev=Vessels{1,i-1}{1,p_i}(1,:)';
        ul=(x0-x0_prev)/norm(x0-x0_prev);
        
        % Determine normal vector to to plane defined by previous 3 points
        if i==2
            un=[1; 0; 0];
        else
            if rem(p_i,2)==1 % odd
                a_i=p_i+1;      % "Aunt" index
                gp_i=(p_i+1)/2; % "Grandparent" index
            else % even
                a_i=p_i-1;	% "Aunt" index
                gp_i=p_i/2;	% "Grandparent" index
            end
            A=Vessels{1,i-2}{1,gp_i}(1,:)';
            B=Vessels{1,i-1}{1,a_i}(2,:)';
            C=x0;

            AB=B-A;
            AC=C-A;
            N=cross(AB,AC);
            un=N/norm(N);
        end
        
        % Find x_tilde point, before rotation by theta
        R_psi = rotation_matrix(un,psi);
        x_tilde=R_psi*current_length*ul+x0;
        
        % Compute and record endpoint
        R_theta = rotation_matrix(ul,theta);
        x_endpoint=R_theta*(x_tilde-x0)+x0;
        
        Vessels{1,i}{1,j}(2,:)=x_endpoint;
    end
    
end
toc

plot_3D_tubes

end