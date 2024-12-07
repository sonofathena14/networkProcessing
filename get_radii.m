function vessel_radii = get_radii(details, cpt_info, TaperID, Scale, ploton)

global data;

Scale = Scale/10; % convert to cm
num_vessels = height(details);

vessel_radii = cell(num_vessels, 7);
vessel_radii{1, 1} = "VESSEL";
vessel_radii{1, 2} = "LENGTH";
vessel_radii{1, 3} = "RIN";
vessel_radii{1, 4} = "ROUT";
vessel_radii{1, 5} = "ERROR";
vessel_radii{1, 6} = "AREA CONSIDERED";
vessel_radii{1, 7} = "TAPERING COEFFICIENTS";

for i = 2:num_vessels
    vessel_radii{i, 1} = details{i, 1};  % Vessel ID
    vessel_radii{i, 2} = details{i, 3};  % Vessel length

    VesselName = ['Vessel:' ' ' details{i,1}(1,:)];

    radius = 0;
    error  = 0;

    radius_info = details{i,2}(2:end,4);
    length_info = details{i,2}(2:end,5);

    if ploton ==1
        figure(i)
        if isempty(cpt_info{i, 2})
            plot(length_info*Scale,radius_info*Scale,'b*','linewidth',3,'markersize',7);
        else
            plot(length_info*Scale,radius_info*Scale,'b*',cpt_info{i,2}(:,2)*Scale,cpt_info{i,2}(:,1)*Scale,'mo','linewidth',3,'markersize',7);
        end
        title(VesselName,'fontsize',16);
        set(gca,'fontsize',16)
        xlabel('Length (cm)');
        ylabel('Radius (cm)');
    end;
    if isempty(cpt_info{i, 2})
        radii_between  = radius_info;
        length_between = length_info;

        rin   = mean(radius_info);
        rout  = rin;
        error = std(radius_info);
    elseif height(cpt_info{i, 2}) == 1
        cpt_location = cpt_info{i, 2}(1, 2);
        cpt_index = find(length_info == cpt_location);

        radii_between  = radius_info(cpt_index:end);
        length_between = length_info(cpt_index:end);

        rin  = mean(radii_between);
        rout = rin;
        
        error  = std(radii_between);   
    else % more than one change-point
        slopes_vector     = cpt_info{i, 3};
        changepoint_radii = cpt_info{i, 2}(:, 2);
        changepoint_radii = sort(changepoint_radii);

        if ismember(details{i, 1},TaperID)  % Vessel will taper
            for j = 1:length(cpt_info{i,2})
                ID(j) = find(length_info == cpt_info{i,2}(j,2));
            end

            length_between = length_info(ID(1):ID(end));
            radii_between  = radius_info(ID(1):ID(end));

            % Using that information, find which changepoints that slope is in
            % between
            rin_est   = radii_between(1);
            rout_est  = radii_between(end);

            data.x = (length_between-length_between(1))*Scale;
            data.r = radii_between*Scale;

            k1 = (rin_est-rout_est)*Scale;
            k2 = 1;
            k3 = rout_est*Scale;
            pars = [k1 k2 k3];

            rad = k1*exp(-k2*data.x)+k3;
            opts = optimset('MaxIter',5000,'MaxFunEvals',5000);
            [xopt, ~, ~, ~] = fminsearch(@model_fmin,pars,opts);
            k1opt(i)   = xopt(1);
            k2opt(i)   = xopt(2);
            k3opt(i)   = xopt(3);
            rad        = k1opt(i)*exp(-k2opt(i)*data.x)+k3opt(i);
            rad_ext    = k1opt(i)*exp(-k2opt(i)*length_info*Scale)+k3opt(i);
            rin        = rad(1)/Scale;
            rout       = rad(end)/Scale;
            error      = std(radii_between);
            N = length(radius_info)-length(radii_between);
            rad_mod = [radii_between' radii_between(end)*ones(1,N)]';
        else
            % Find the min of the abs value of slopes vector and its index
            % after first changepoint
            flattest_slope = min(abs(slopes_vector(2:end)));
            flattest_index = find(abs(slopes_vector) == flattest_slope);
            if length(flattest_index) > 1
                %fprintf("Length of flattest index is > 1 at %d\n", i)
                flattest_index = flattest_index(2);
            end

            % Using that information, find which changepoints that slope is in
            % between
            index_before = 0;
            index_after  = 0;
            if flattest_index == height(slopes_vector)
                index_before = find(length_info == changepoint_radii(end));
                index_after  = height(length_info);
            else
                index_before = find(length_info == changepoint_radii(flattest_index-1));
                index_after  = find(length_info == changepoint_radii(flattest_index));
            end

            % Find the indices in the radii vector in details where
            % those changepoints are
            % Determine what percentage of the data is in between those two
            % points
            percent_between = (index_after-index_before+1)/height(length_info);

            % If the segment covers > 25% of the vessel data
            if percent_between >= .25
                % Make a new vector consisting of only the data points captured
                % between the two changepoints
                radii_between  = radius_info(index_before:index_after);
                length_between = length_info(index_before:index_after);

                rin  = mean(radii_between); % radius = avg of this vector
                rout = rin;
                error= std(radii_between);% error = std of this vector
            else
                % Find the index of the first changepoint
                first_cpt       = changepoint_radii(1);
                first_cpt_index = find(length_info == first_cpt);

                % Determine what percentage of the data is from the beginning
                % to the first changepoint
                percent_between = first_cpt_index/height(length_info);

                % If this percentage is > 25% of the vessel data
                if percent_between >= .25
                    % If last two elements of slope vector are far enough from
                    % each other (threshold to be determined)
                    last_two_slopes_diff = (slopes_vector(end)-slopes_vector(end-1))/slopes_vector(end-1);
                    if abs(last_two_slopes_diff) >= 2.8
                        last_cpt       = changepoint_radii(end);
                        last_cpt_index = find(length_info == last_cpt);

                        % Create a vector consisting of the data points between
                        % the first changepoint and the last changepoint
                        radii_between  = radius_info(first_cpt_index:last_cpt_index);
                        length_between = length_info(first_cpt_index:last_cpt_index);

                        rin   = mean(radii_between); % radius = avg of this vector
                        rout  = rin;
                        error = std(radii_between); % error = std of this vector
                    else
                        % Create a vector consisting the data points between
                        % the first changepoint and the end of radius vector
                        radii_between  = radius_info(first_cpt_index:end);
                        length_between = length_info(first_cpt_index:end);

                        rin   = mean(radii_between);  % radius = avg of this vector
                        rout  = rin;
                        error = std(radii_between);   % error = std of this vector
                    end
                else
                    % Find the index of the second changepoint
                    second_cpt       = changepoint_radii(2);
                    second_cpt_index = find(length_info == second_cpt);

                    % Determine what percentage of the data is from the
                    % beginning to the second changepoint
                    percent_between = second_cpt_index/height(length_info);

                    % If this percentage is > 50% of the vessel data
                    if percent_between >= .50
                        % Find the index closest to the 25% point of the
                        % data above
                        pct_25_index = ceil(.25 * height(length_info));
                        % If last two elements of slope vector are far enough from
                        % each other (threshold to be determined)
                        last_two_slopes_diff = (slopes_vector(end)-slopes_vector(end-1))/slopes_vector(end-1);
                        if abs(last_two_slopes_diff) >= 2.8
                            last_cpt       = changepoint_radii(end);
                            last_cpt_index = find(length_info == last_cpt);

                            % Create a vector consisting of the data points between
                            % the 25% index and the last changepoint
                            radii_between  = radius_info(pct_25_index:last_cpt_index);
                            length_between = length_info(pct_25_index:last_cpt_index);

                            rin   = mean(radii_between); % radius = avg of this vector
                            rout  = rin;
                            error = std(radii_between); % error = std of this vector
                        else
                            % Create a vector consisting the data points between
                            % the 25% index and the end of radius vector
                            radii_between  = radius_info(pct_25_index:end);
                            length_between = length_info(pct_25_index:end);

                            rin   = mean(radii_between); % radius = avg of this vector
                            rout  = rin;
                            error = std(radii_between); % error = std of this vector
                        end
                    else
                        % If last two elements of slope vector are far enough from
                        % each other (threshold to be determined)
                        last_two_slopes_diff = (slopes_vector(end)-slopes_vector(end-1))/slopes_vector(end-1);
                        if abs(last_two_slopes_diff) >= 2.8
                            last_cpt       = changepoint_radii(end);
                            last_cpt_index = find(length_info == last_cpt);

                            % Create a vector consisting of the data points between
                            % the second changepoint and the last changepoint
                            radii_between  = radius_info(second_cpt_index:last_cpt_index);
                            length_between = length_info(second_cpt_index:last_cpt_index);

                            rin   = mean(radii_between); % radius = avg of this vector
                            rout  = rin;
                            error = std(radii_between); % error = std of this vector
                        else
                            % Create a vector consisting the data points between
                            % the second changepoint and the end of radius vector
                            radii_between  = radius_info(second_cpt_index:end);
                            length_between = length_info(second_cpt_index:end);

                            rin   = mean(radii_between); % radius = avg of this vector
                            rout  = rin;
                            error = std(radii_between); % error = std of this vector
                        end
                    end
                end
            end
        end
    end
    if ploton == 1
        if ismember(details{i, 1},TaperID)
            figure(i); hold on;
            plot(length_between*Scale,rad,'m-','linewidth',3,'markersize',7);
            plot(length_between*Scale,rad+error*Scale*ones(size(length_between)),'m--','linewidth',3,'markersize',7);
            plot(length_between*Scale,rad-error*Scale*ones(size(length_between)),'m--','linewidth',3,'markersize',7);
            set(gca,'fontsize',16);
            xlabel('length (cm)');
            ylabel('radius (cm)');
            legend('Vessel radii', 'Est seg radii','Est vessel radii','Tapered radii')

            figure(i+100); hold on;
            plot(length_info*Scale,radius_info*Scale,'b*',cpt_info{i,2}(:,2)*Scale,cpt_info{i,2}(:,1)*Scale,'mo','linewidth',3,'markersize',7);
            plot(length_info*Scale,rad_mod*Scale,'c*',length_info*Scale,rad_ext,'m-','linewidth',3,'markersize',7);
            plot(length_info*Scale,rad_ext+error*Scale*ones(size(length_info)),'m--','linewidth',3,'markersize',7);
            plot(length_info*Scale,rad_ext-error*Scale*ones(size(length_info)),'m--','linewidth',3,'markersize',7);
            set(gca,'fontsize',16);
            xlabel('length (cm)');
            ylabel('radius (cm)');
            legend('Vessel radii','Est vessel radii','Tapered radii')
            xlim([0 length_info(end)*Scale]);
            pause;
        else
            figure(i);hold on;
            plot(length_info*Scale,rin*Scale*ones(size(length_info)),'m','LineWidth',3)
            plot(length_info*Scale,(rin+error)*Scale*ones(size(length_info)),'m--','LineWidth',3)
            plot(length_info*Scale,(rin-error)*Scale*ones(size(length_info)),'m--','LineWidth',3)
            xlim([0 length_info(end)*Scale]);
        end
        pause;
    end
    vessel_radii{i, 2} = vessel_radii{i, 2}*Scale;
    vessel_radii{i, 3} = rin*Scale;
    vessel_radii{i, 4} = rout*Scale;
    vessel_radii{i, 5} = error*Scale;
    vessel_radii{i, 6} = [length_between*Scale radii_between*Scale];
    if ismember(details{i, 1},TaperID)
        ID = str2num(details{i, 1});
        vessel_radii{i, 7} = [ID k1opt(i) k2opt(i) k3opt(i)];
    end
end
end % function %
       