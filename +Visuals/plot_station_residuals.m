function plot_station_residuals(Measurement_Table, station_names, P_zz)
    % Measurement_Table: The table containing residuals and station_id
    % station_names: Cell array of strings, e.g., {'Atoll', 'Diego Garcia', ...}

    % Calculate Residuals
    res_range = Measurement_Table.apparent_range_meters - Measurement_Table.computed_range_meters;
    res_rate  = Measurement_Table.apparent_range_rate_meters_s - Measurement_Table.computed_range_rate_meters_s;
    t_sec     = Measurement_Table.time_sec_past_epoch;
    u_stations = unique(Measurement_Table.station_id);

    three_sigma_range_meas = squeeze(P_zz(1,1,:)) * 3;
    three_sigma_range_rate_meas = squeeze(P_zz(2,2,:)) * 3;
    
    % Setup Figure
    figure('Color', 'w');
    colors = lines(length(u_stations)); % Distinct color palette
    
    % --- Top Subplot: Range Residuals ---
    hold on; grid on;
    for k = 1:length(u_stations)
        idx = (Measurement_Table.station_id == u_stations(k));
        plot(t_sec(idx), res_range(idx), '.', 'MarkerFaceColor', colors(k,:), ...
            'Color', colors(k,:), 'DisplayName', station_names{u_stations(k)});
    end

    hold on;

    plot(t_sec, three_sigma_range_meas,'k-');
    plot(t_sec, -1 * three_sigma_range_meas,'k-')

    ylim([-200,200]);
    ylabel('Range Residual (meters)'); title('Station-Dependent Range Residuals');
    legend('Location', 'bestoutside'); set(gca, 'FontSize', 12);

    % --- Bottom Subplot: Range-Rate Residuals ---
    figure('Color', 'w');
    hold on; grid on;
    for k = 1:length(u_stations)
        idx = (Measurement_Table.station_id == u_stations(k));
        plot(t_sec(idx), res_rate(idx), '.', 'MarkerFaceColor', colors(k,:), ...
            'Color', colors(k,:), 'DisplayName', station_names{u_stations(k)});
    end
    
    hold on;

    plot(t_sec, three_sigma_range_rate_meas,'k-');
    plot(t_sec, -1 * three_sigma_range_rate_meas,'k-')

    ylim([-1,1]);
    ylabel('Range-Rate Residual (meters/s)'); xlabel('Time since epoch (s)');
    title('Station-Dependent Range-Rate Residuals');
    legend('Location', 'bestoutside'); set(gca, 'FontSize', 12);
end