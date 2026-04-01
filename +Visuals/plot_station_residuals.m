function plot_station_residuals(Measurement_Table, station_names)
    % Measurement_Table: The table containing residuals and station_id
    % station_names: Cell array of strings, e.g., {'Atoll', 'Diego Garcia', ...}

    % Calculate Residuals
    res_range = Measurement_Table.apparent_range_km - Measurement_Table.computed_range_km;
    res_rate  = Measurement_Table.apparent_range_rate_km_s - Measurement_Table.computed_range_rate_km_s;
    t_sec     = Measurement_Table.time_sec_past_epoch;
    u_stations = unique(Measurement_Table.station_id);
    
    % Setup Figure
    figure('Color', 'w', 'Position', [100, 100, 1000, 800]);
    colors = lines(length(u_stations)); % Distinct color palette
    
    % --- Top Subplot: Range Residuals ---
    subplot(2,1,1); hold on; grid on;
    for k = 1:length(u_stations)
        idx = (Measurement_Table.station_id == u_stations(k));
        plot(t_sec(idx), res_range(idx), 'o', 'MarkerFaceColor', colors(k,:), ...
            'Color', colors(k,:), 'DisplayName', station_names{u_stations(k)});
    end
    ylabel('Range Residual (km)'); title('Station-Dependent Range Residuals');
    legend('Location', 'bestoutside'); set(gca, 'FontSize', 12);

    % --- Bottom Subplot: Range-Rate Residuals ---
    subplot(2,1,2); hold on; grid on;
    for k = 1:length(u_stations)
        idx = (Measurement_Table.station_id == u_stations(k));
        plot(t_sec(idx), res_rate(idx), 's', 'MarkerFaceColor', colors(k,:), ...
            'Color', colors(k,:), 'DisplayName', station_names{u_stations(k)});
    end
    ylabel('Range-Rate Residual (km/s)'); xlabel('Time since epoch (s)');
    title('Station-Dependent Range-Rate Residuals');
    legend('Location', 'bestoutside'); set(gca, 'FontSize', 12);
end