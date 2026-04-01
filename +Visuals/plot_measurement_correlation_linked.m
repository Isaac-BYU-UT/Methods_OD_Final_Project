function plot_measurement_correlation_linked(Measurement_Table)
    % Plots Truth vs Model as discrete points with vertical connector lines
    t_sec = Measurement_Table.time_sec_past_epoch;
    obs_r = Measurement_Table.apparent_range_km;
    comp_r = Measurement_Table.computed_range_km;
    obs_rr = Measurement_Table.apparent_range_rate_km_s;
    comp_rr = Measurement_Table.computed_range_rate_km_s;
    
    figure('Color', 'w', 'Name', 'Measurement Correlation - Linked');
    
    % --- 1. Range Correlation ---
    subplot(2,1,1); hold on; grid on;
    
    % Draw Vertical "Residual" Lines first so they sit behind the points
    % This uses a single plot call with NaNs to be computationally efficient
    line_x = [t_sec'; t_sec'; NaN(1, length(t_sec))];
    line_y_r = [obs_r'; comp_r'; NaN(1, length(t_sec))];
    plot(line_x(:), line_y_r(:), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    % Plot Discrete Points
    plot(t_sec, obs_r, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Observed (o)');
    plot(t_sec, comp_r, 'rd', 'MarkerSize', 6, 'DisplayName', 'Computed (d)');
    
    ylabel('Range (km)'); 
    title('Range: Observed vs Computed (Lines represent Residuals)');
    legend('show', 'Location', 'best');

    % --- 2. Range-Rate Correlation ---
    subplot(2,1,2); hold on; grid on;
    
    % Draw Vertical "Residual" Lines
    line_y_rr = [obs_rr'; comp_rr'; NaN(1, length(t_sec))];
    plot(line_x(:), line_y_rr(:), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    % Plot Discrete Points
    plot(t_sec, obs_rr, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Observed (o)');
    plot(t_sec, comp_rr, 'bd', 'MarkerSize', 6, 'DisplayName', 'Computed (d)');
    
    ylabel('Range-Rate (km/s)'); xlabel('Time (s)');
    title('Doppler: Observed vs Computed');
    legend('show', 'Location', 'best');
end