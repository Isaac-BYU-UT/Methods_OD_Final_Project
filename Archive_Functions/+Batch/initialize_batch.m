function [S, batch] = initialize_batch(S, ENV, batch)
% Initialize batch least squares OD structure
% Following similar conventions as EKF initialization

    batch.N_states = 6; % Position and velocity (exclude drag for now)
    batch.station_meas_mask = [S.Scenario.Atoll_on, S.Scenario.Diego_Garcia_on, S.Scenario.Arecibo_on];

    % Set measurement covariances based on scenario settings
    if (~S.Scenario.range_on)
        S.Stations(1).Covariance(1,1) = 1e16;
        S.Stations(2).Covariance(1,1) = 1e16;
        S.Stations(3).Covariance(1,1) = 1e16;
    end

    if (~S.Scenario.range_rate_on)
        S.Stations(1).Covariance(2,2) = 1e16;
        S.Stations(2).Covariance(2,2) = 1e16;
        S.Stations(3).Covariance(2,2) = 1e16;
    end

    % Initial state
    X0_nominal = [ ...
        S.IC_Sat_Epoch.position_ECI_meters;
        S.IC_Sat_Epoch.velocity_ECI_meters_per_second;
    ];

    batch.X0_nominal = X0_nominal;
    batch.P_cov_0 = S.StateCovariances.P_Covariance_States(1:6, 1:6); % a priori covariance
    
    % Observation times and count
    batch.t_obs = S.ref_data.Actual_Measurements.time_sec_past_epoch;
    batch.N_obs = length(batch.t_obs);

    % Pre-allocate storage for measurements and residuals
    batch.Y_obs = zeros(2, batch.N_obs);          % Observed measurements
    batch.Y_computed = zeros(2, batch.N_obs);     % Computed measurements
    batch.Y_residuals = zeros(2, batch.N_obs);    % Residuals (O-C)

    % Propagation settings
    batch.X_input = X0_nominal;

end
