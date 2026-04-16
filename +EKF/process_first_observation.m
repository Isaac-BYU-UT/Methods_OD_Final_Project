function ekf = process_first_observation(ekf, S, ENV)

    i = ekf.current_index; % THIS SHOULD ALWAYS BE 1!
    X_nominal = ekf.X_input(1:ekf.N_states); % This is just our initial condition
    [curr_meas, ekf] = EKF.compute_measurement(ekf, S, ENV, X_nominal); % Before updating X

    [X_updated, P_cov_updated, dx] = EKF.ekf_update( ...
        ekf, X_nominal, ekf.P_cov, curr_meas, S, ENV);
    
    r_sat_updated_ECI_m = X_updated(1:3);
    v_sat_updated_ECI_m_s = X_updated(4:6);

    ekf.Y_postfit(:,i) = Measurements.Compute_Range_Range_Rate( ...
        r_sat_updated_ECI_m, v_sat_updated_ECI_m_s, curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s);
    
    EKF.log_step(ekf,dx,P_cov_updated,i);   % LOG before we overwrite ekf

    STM_reset = eye(ekf.N_states);
    ekf.X_input = [X_updated; STM_reset(:)]; % Prep for next iteration
    ekf.P_cov   = P_cov_updated;             % Prep for next iteration

    
end