function ekf = process_first_observation(ekf, S, ENV)

    i = 1;

    meas = compute_measurement(ekf, S, ENV, ekf.X_input(1:ekf.N_states), i);

    [X_updated, P_cov, ~] = ekf_update( ...
        ekf, ekf.X_input(1:ekf.N_states), ekf.P_cov, meas, S, ENV);

    ekf.Y_postfit(:,1) = compute_observation( ...
        X_updated, meas.r_stn_ECI_m, meas.v_stn_ECI_m_s);

    STM_reset = zeros(ekf.N_states, 1);
    ekf.X_input = [X_updated; STM_reset(:)];
    ekf.P_cov   = P_cov;
end