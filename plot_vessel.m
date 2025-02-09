function output = plot_vessel(Name, Scale, pv)
close all;

vessel_detail = load(strcat('FluidsInput/Changepoint_',Name,'.mat'), 'vessel_details');
details = vessel_detail.vessel_details;
changepoints = load(strcat('FluidsInput/Changepoint_',Name,'.mat'), 'changepoint_location');
cpt_info = changepoints.changepoint_location;
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

all_ids = details(:,1); % Convert first column to numeric array
all_ids_num = cellfun(@str2double, all_ids);
ves_id = find(all_ids_num == pv);
daughters = details{ves_id, 8};
daughter_id = daughters{1};
vessels = [pv daughter_id];

for i = 1:length(vessels)
    all_ids = details(:,1); % Convert first column to numeric array
    all_ids_num = cellfun(@str2double, all_ids);
    ves_id = find(all_ids_num == vessels(i));

    vessel_radii{ves_id, 1} = details{ves_id, 1};  % Vessel ID
    vessel_radii{ves_id, 2} = details{ves_id, 3};  % Vessel length

    VesselName = ['Vessel:' ' ' details{ves_id,1}(1,:)];

    radius = 0;
    error  = 0;

    radius_info = details{ves_id,2}(2:end,4);
    length_info = details{ves_id,2}(2:end,5);

    figure(i)
    if isempty(cpt_info{ves_id, 2})
        plot(length_info*Scale,radius_info*Scale,'b*','linewidth',3,'markersize',7);
    else
        plot(length_info*Scale,radius_info*Scale,'b*',cpt_info{ves_id,2}(:,2)*Scale,cpt_info{ves_id,2}(:,1)*Scale,'mo','linewidth',3,'markersize',7);
    end
    title(VesselName,'fontsize',16);
    set(gca,'fontsize',16)
    xlabel('Length (cm)');
    ylabel('Radius (cm)');
    if isempty(cpt_info{ves_id, 2})
        radii_between  = radius_info;
        length_between = length_info;

        rin   = mean(radius_info);
        rout  = rin;
        error = std(radius_info);
    elseif height(cpt_info{ves_id, 2}) == 1
        cpt_location = cpt_info{ves_id, 2}(1, 2);
        cpt_index = find(length_info == cpt_location);

        radii_between  = radius_info(cpt_index:end);
        length_between = length_info(cpt_index:end);

        rin  = mean(radii_between);
        rout = rin;
        
        error  = std(radii_between);   
    else % more than one change-point
        slopes_vector     = cpt_info{ves_id, 3};
        changepoint_radii = cpt_info{ves_id, 2}(:, 2);
        changepoint_radii = sort(changepoint_radii);
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
    figure(i);hold on;
    plot(length_info*Scale,rin*Scale*ones(size(length_info)),'m','LineWidth',3)
    plot(length_info*Scale,(rin+error)*Scale*ones(size(length_info)),'m--','LineWidth',3)
    plot(length_info*Scale,(rin-error)*Scale*ones(size(length_info)),'m--','LineWidth',3)
    xlim([0 length_info(end)*Scale]);
end
output = 1;
end