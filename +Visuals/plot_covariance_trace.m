function plot_covariance_trace(ekf)
    % Plots the trace of the state covariance matrix over time
    % Shows how uncertainty evolves after each propagation and measurement update
    %
    % Inputs:
    %   ekf - EKF structure containing trace history arrays
    
    if ~isfield(ekf, 'trace_post_propagation') || ~isfield(ekf, 'trace_post_update')
        error('EKF structure does not contain covariance trace data. Ensure log_step was called during filter execution.');
    end
    
    % Convert time to datetime or just use index
    t_obs_hours = ekf.t_obs / 3600;  % Convert seconds to hours
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
    
    % --- Plot 1: Trace over time with propagation and update steps ---
    subplot(1, 2, 1); hold on;
    
    % Remove zero entries (unfilled slots from initialization)
    valid_idx = ekf.t_obs > 0;
    
    plot(t_obs_hours(valid_idx), ekf.trace_post_propagation(valid_idx), 'b-o', ...
        'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'After Propagation (P_bar)');
    plot(t_obs_hours(valid_idx), ekf.trace_post_update(valid_idx), 'r-s', ...
        'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'After Update (P)');
    plot(t_obs_hours(valid_idx), ekf.trace_pre_propagation(valid_idx), 'k-^', ...
        'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'Before Propagation', 'LineStyle', '--');
    
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Covariance Trace (m^2 + (m/s)^2)', 'FontSize', 12);
    title('State Covariance Trace Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    
    % --- Plot 2: Change in trace through each step ---
    subplot(1, 2, 2); hold on;
    
    % Calculate differences
    prop_change = ekf.trace_post_propagation(valid_idx) - ekf.trace_pre_propagation(valid_idx);
    update_change = ekf.trace_post_update(valid_idx) - ekf.trace_post_propagation(valid_idx);
    
    % Create bar chart
    x = 1:sum(valid_idx);
    width = 0.35;
    
    bar(x - width/2, prop_change, width, 'FaceColor', [0.2 0.4 0.8], ...
        'DisplayName', 'Propagation Step (Increase)', 'EdgeColor', 'k', 'LineWidth', 0.5);
    bar(x + width/2, update_change, width, 'FaceColor', [0.8 0.2 0.2], ...
        'DisplayName', 'Update Step (Decrease)', 'EdgeColor', 'k', 'LineWidth', 0.5);
    
    xlabel('Observation Index', 'FontSize', 12);
    ylabel('Change in Trace (m^2 + (m/s)^2)', 'FontSize', 12);
    title('Trace Change Per Filter Step', 'FontSize', 13, 'FontWeight', 'bold');
    yline(0, 'k--', 'LineWidth', 1);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    
    % Summary statistics
    fprintf('\n--- Covariance Trace Statistics ---\n');
    fprintf('Initial trace (before propagation):  %.6e\n', ekf.trace_pre_propagation(find(valid_idx, 1)));
    fprintf('Final trace (after last update):     %.6e\n', ekf.trace_post_update(find(valid_idx, 1, 'last')));
    fprintf('Minimum trace during filter:        %.6e\n', min(ekf.trace_post_update(valid_idx)));
    fprintf('Maximum trace during filter:        %.6e\n', max(ekf.trace_post_propagation(valid_idx)));
    fprintf('Mean propagation step change:       %.6e\n', mean(prop_change));
    fprintf('Mean update step change:            %.6e\n', mean(update_change));
    fprintf('Total uncertainty reduction:        %.2f%%\n', ...
        100 * (ekf.trace_pre_propagation(find(valid_idx, 1)) - ekf.trace_post_update(find(valid_idx, 1, 'last'))) / ...
        ekf.trace_pre_propagation(find(valid_idx, 1)));
end
