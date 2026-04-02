function plot_station_residuals_3_sigma(Measurement_Table, station_names, IC)
    % Calculate Residuals
    res_range = Measurement_Table.apparent_range_km - Measurement_Table.computed_range_km;
    res_rate  = Measurement_Table.apparent_range_rate_km_s - Measurement_Table.computed_range_rate_km_s;
    t_sec     = Measurement_Table.time_sec_past_epoch;
    u_stations = unique(Measurement_Table.station_id);
    
    figure('Color', 'w');
    colors = lines(length(u_stations)); 
    
    % --- Top Subplot: Range Residuals ---
    subplot(2,1,1); hold on; grid on;
    for k = 1:length(u_stations)
        stn_id = u_stations(k);
        idx = (Measurement_Table.station_id == stn_id);
        
        % Plot Residuals
        plot(t_sec(idx), res_range(idx), 'o', 'MarkerFaceColor', colors(k,:), ...
            'Color', colors(k,:), 'DisplayName', station_names{stn_id});
        
        % Calculate 3-sigma (sqrt of variance * 3)
        sigma_range = sqrt(IC.Stations(stn_id).Covariance(1,1));
        three_sigma = 3 * sigma_range;
        
        % Plot Bounds (Plotting as a line across the data range)
        line([min(t_sec(idx)) max(t_sec(idx))], [three_sigma three_sigma], ...
            'Color', colors(k,:), 'LineStyle', '--', 'HandleVisibility', 'off');
        line([min(t_sec(idx)) max(t_sec(idx))], [-three_sigma -three_sigma], ...
            'Color', colors(k,:), 'LineStyle', '--', 'HandleVisibility', 'off');
    end
    ylabel('Range Residual (km)'); title('Range Residuals with 3\sigma Bounds');
    legend('Location', 'bestoutside'); set(gca, 'FontSize', 12);

    % --- Bottom Subplot: Range-Rate Residuals ---
    subplot(2,1,2); hold on; grid on;
    for k = 1:length(u_stations)
        stn_id = u_stations(k);
        idx = (Measurement_Table.station_id == stn_id);
        
        plot(t_sec(idx), res_rate(idx), 's', 'MarkerFaceColor', colors(k,:), ...
            'Color', colors(k,:), 'DisplayName', station_names{stn_id});

        % Calculate 3-sigma
        sigma_rate = sqrt(IC.Stations(stn_id).Covariance(2,2));
        three_sigma_rate = 3 * sigma_rate;

        line([min(t_sec(idx)) max(t_sec(idx))], [three_sigma_rate three_sigma_rate], ...
            'Color', colors(k,:), 'LineStyle', '--', 'HandleVisibility', 'off');
        line([min(t_sec(idx)) max(t_sec(idx))], [-three_sigma_rate -three_sigma_rate], ...
            'Color', colors(k,:), 'LineStyle', '--', 'HandleVisibility', 'off');
    end
    ylabel('Range-Rate Residual (km/s)'); xlabel('Time since epoch (s)');
    title('Range-Rate Residuals with 3\sigma Bounds');
    legend('Location', 'bestoutside'); set(gca, 'FontSize', 12);
end