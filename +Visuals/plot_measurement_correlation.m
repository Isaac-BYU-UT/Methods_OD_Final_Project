function plot_measurement_correlation(Measurement_Table)
    % Plots Truth (Observed) vs Model (Computed) to find biases
    t_sec = Measurement_Table.time_sec_past_epoch;
    
    figure('Color', 'w', 'Name', 'Measurement Correlation');
    
    % Range Correlation
    subplot(2,1,1); hold on; grid on;
    plot(t_sec, Measurement_Table.apparent_range_km, 'k.', 'MarkerSize', 12, 'DisplayName', 'Observed (Truth)');
    plot(t_sec, Measurement_Table.computed_range_km, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Computed (Model)');
    ylabel('Range (km)'); title('Range Correlation: Observed vs Computed');
    legend('show');

    % Range-Rate Correlation
    subplot(2,1,2); hold on; grid on;
    plot(t_sec, Measurement_Table.apparent_range_rate_km_s, 'k.', 'MarkerSize', 12, 'DisplayName', 'Observed (Truth)');
    plot(t_sec, Measurement_Table.computed_range_rate_km_s, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Computed (Model)');
    ylabel('Range-Rate (km/s)'); xlabel('Time (s)');
    title('Doppler Correlation: Observed vs Computed');
    legend('show');
end