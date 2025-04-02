function output = plot_vessel2(Name, Scale, pv)
close all;

vessel_detail = load(strcat('Output/Vessels_',Name,'.mat'), 'new_vessel_details');
details = vessel_detail.new_vessel_details;
changepoints = load(strcat('Output/Vessels_',Name,'.mat'), 'changepoint_location');
cpt_info = changepoints.changepoint_location;

all_ids = details(:,1); % Convert first column to numeric array
all_ids_num = cellfun(@str2double, all_ids);
ves_id = find(all_ids_num == pv);
daughters = details{ves_id, 8};
daughter_id = daughters{1};
vessels = [pv daughter_id.'];

% Create a larger figure with a grid layout for the subplots
figure;

% Loop over vessels and plot each in a separate subplot
for i = 1:length(vessels)
    all_ids = details(:,1); % Convert first column to numeric array
    all_ids_num = cellfun(@str2double, all_ids);
    ves_id = find(all_ids_num == vessels(i));
        
    vessel_radii{ves_id, 1} = details{ves_id, 1};  % Vessel ID
    vessel_radii{ves_id, 2} = details{ves_id, 3};  % Vessel length
    names = ["Parent" "Daughter 1" "Daughter 2"];

    VesselName = strcat('Vessel:', ' ', details{ves_id, 1}(1, :), ' ', names(:, i));

    radius = 0;
    error  = 0;

    radius_info = details{ves_id,2}(2:end,4);
    length_info = details{ves_id,2}(2:end,5);

    % Set the current subplot
    subplot(1, 3, mod(i-1, 3) + 1); % Adjust number of subplots if needed (1 row, 3 columns)

    if isempty(cpt_info{ves_id, 2})
        plot(length_info*Scale, radius_info*Scale, 'b*', 'LineWidth', 3, 'MarkerSize', 7);
    else
        plot(length_info*Scale, radius_info*Scale, 'b*', cpt_info{ves_id, 2}(:,2)*Scale, cpt_info{ves_id, 2}(:,1)*Scale, 'mo', 'LineWidth', 3, 'MarkerSize', 7);
    end
    title(VesselName, 'FontSize', 16);
    set(gca, 'FontSize', 16);
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
        flattest_slope = min(abs(slopes_vector(2:end)));
        flattest_index = find(abs(slopes_vector) == flattest_slope);
        if length(flattest_index) > 1
            flattest_index = flattest_index(2);
        end

        index_before = 0;
        index_after  = 0;
        if flattest_index == height(slopes_vector)
            index_before = find(length_info == changepoint_radii(end));
            index_after  = height(length_info);
        else
            index_before = find(length_info == changepoint_radii(flattest_index-1));
            index_after  = find(length_info == changepoint_radii(flattest_index));
        end

        percent_between = (index_after - index_before + 1) / height(length_info);

        if percent_between >= .25
            radii_between  = radius_info(index_before:index_after);
            length_between = length_info(index_before:index_after);

            rin  = mean(radii_between);
            rout = rin;
            error= std(radii_between);
        else
            first_cpt       = changepoint_radii(1);
            first_cpt_index = find(length_info == first_cpt);

            percent_between = first_cpt_index / height(length_info);

            if percent_between >= .25
                last_two_slopes_diff = (slopes_vector(end) - slopes_vector(end - 1)) / slopes_vector(end - 1);
                if abs(last_two_slopes_diff) >= 2.8
                    last_cpt       = changepoint_radii(end);
                    last_cpt_index = find(length_info == last_cpt);

                    radii_between  = radius_info(first_cpt_index:last_cpt_index);
                    length_between = length_info(first_cpt_index:last_cpt_index);

                    rin   = mean(radii_between);
                    rout  = rin;
                    error = std(radii_between);
                else
                    radii_between  = radius_info(first_cpt_index:end);
                    length_between = length_info(first_cpt_index:end);

                    rin   = mean(radii_between);
                    rout  = rin;
                    error = std(radii_between);
                end
            else
                second_cpt       = changepoint_radii(2);
                second_cpt_index = find(length_info == second_cpt);

                percent_between = second_cpt_index / height(length_info);

                if percent_between >= .50
                    pct_25_index = ceil(.25 * height(length_info));

                    last_two_slopes_diff = (slopes_vector(end) - slopes_vector(end - 1)) / slopes_vector(end - 1);
                    if abs(last_two_slopes_diff) >= 2.8
                        last_cpt       = changepoint_radii(end);
                        last_cpt_index = find(length_info == last_cpt);

                        radii_between  = radius_info(pct_25_index:last_cpt_index);
                        length_between = length_info(pct_25_index:last_cpt_index);

                        rin   = mean(radii_between);
                        rout  = rin;
                        error = std(radii_between);
                    else
                        radii_between  = radius_info(pct_25_index:end);
                        length_between = length_info(pct_25_index:end);

                        rin   = mean(radii_between);
                        rout  = rin;
                        error = std(radii_between);
                    end
                else
                    last_two_slopes_diff = (slopes_vector(end) - slopes_vector(end - 1)) / slopes_vector(end - 1);
                    if abs(last_two_slopes_diff) >= 2.8
                        last_cpt       = changepoint_radii(end);
                        last_cpt_index = find(length_info == last_cpt);

                        radii_between  = radius_info(second_cpt_index:last_cpt_index);
                        length_between = length_info(second_cpt_index:last_cpt_index);

                        rin   = mean(radii_between);
                        rout  = rin;
                        error = std(radii_between);
                    else
                        radii_between  = radius_info(second_cpt_index:end);
                        length_between = length_info(second_cpt_index:end);

                        rin   = mean(radii_between);
                        rout  = rin;
                        error = std(radii_between);
                    end
                end
            end
        end
    end
    hold on;
    plot(length_info*Scale, rin*Scale*ones(size(length_info)), 'm', 'LineWidth', 3);
    plot(length_info*Scale, (rin + error)*Scale*ones(size(length_info)), 'm--', 'LineWidth', 3);
    plot(length_info*Scale, (rin - error)*Scale*ones(size(length_info)), 'm--', 'LineWidth', 3);
    xlim([0 length_info(end)*Scale]);
    ylim([0 11]);
end

end