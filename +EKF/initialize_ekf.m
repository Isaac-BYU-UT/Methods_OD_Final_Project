function ekf = initialize_ekf(S, ENV)

    ekf.N_states = 7;
    ekf.station_meas_mask = [S.Scenario.Atoll_on, S.Scenario.Diego_Garcia_on, S.Scenario.Arecibo_on];

    if (~S.Scenario.range_on)
        S.Stations(1).Covariance(1,1) = 1e12; % effectively infinite variance for range
        S.Stations(2).Covariance(1,1) = 1e12;
        S.Stations(3).Covariance(1,1) = 1e12;
    end

    if (~S.Scenario.range_rate_on)
        S.Stations(1).Covariance(2,2) = 1e12; % effectively infinite variance for range rate
        S.Stations(2).Covariance(2,2) = 1e12;
        S.Stations(3).Covariance(2,2) = 1e12;
    end

    STM0 = eye(ekf.N_states);
    C_drag_estimate = 1.88;

    X0 = [ ...
        S.IC_Sat_Epoch.position_ECI_meters;
        S.IC_Sat_Epoch.velocity_ECI_meters_per_second;
        C_drag_estimate;
        STM0(:)
    ];

    ekf.time_struct = ENV.time_struct;
    ekf.options = odeset('RelTol',3e-10,'AbsTol',1e-12);

    ekf.P_cov = S.StateCovariances;

    ekf.t_obs = S.ref_data.Actual_Measurements.time_sec_past_epoch;
    ekf.N_obs = length(ekf.t_obs);

    ekf.Y_prefit  = zeros(2, ekf.N_obs);
    ekf.Y_postfit = zeros(2, ekf.N_obs);

    ekf.X_input = X0;
end