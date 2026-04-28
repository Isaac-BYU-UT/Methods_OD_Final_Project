function [S,ekf] = initialize_ekf(S, ENV, ekf)

    ekf.N_states = 6; % Exclude drag
    ekf.station_meas_mask = [S.Scenario.Atoll_on, S.Scenario.Diego_Garcia_on, S.Scenario.Arecibo_on];

    if (~S.Scenario.range_on)
        S.Stations(1).Covariance(1,1) = 1e16; % effectively infinite variance for range
        S.Stations(2).Covariance(1,1) = 1e16;
        S.Stations(3).Covariance(1,1) = 1e16;
    end

    if (~S.Scenario.range_rate_on)
        S.Stations(1).Covariance(2,2) = 1e16; % effectively infinite variance for range rate
        S.Stations(2).Covariance(2,2) = 1e16;
        S.Stations(3).Covariance(2,2) = 1e16;
    end

    STM0 = eye(ekf.N_states);
    % C_drag_estimate = 1.88;

    X0_star = [ ...
        S.IC_Sat_Epoch.position_ECI_meters;
        S.IC_Sat_Epoch.velocity_ECI_meters_per_second;
        % C_drag_estimate;
        STM0(:)
    ];

    ekf.P_cov = S.StateCovariances.P_Covariance_States(1:6,1:6); % a priori p_bar_0, ignore drag
    ekf.t_obs = S.ref_data.Actual_Measurements.time_sec_past_epoch;
    ekf.N_obs = length(ekf.t_obs);

    ekf.Y_prefit  = zeros(2, ekf.N_obs);
    ekf.Y_postfit = zeros(2, ekf.N_obs);
    
    % Covariance trace tracking
    ekf.trace_pre_propagation = zeros(ekf.N_obs, 1);   % trace before propagation
    ekf.trace_post_propagation = zeros(ekf.N_obs, 1);  % trace after propagation (P_bar)
    ekf.trace_post_update = zeros(ekf.N_obs, 1);       % trace after measurement update

    ekf.X_input = X0_star; % this is our X_star_t_minus_1
end