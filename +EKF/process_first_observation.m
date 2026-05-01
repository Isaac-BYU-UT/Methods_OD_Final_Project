function ekf = process_first_observation(ekf, S, ENV)

% Prefit saved in compute_measurement
% Both P_zz and postfit saved in ekf_update

    i = ekf.current_index; % THIS SHOULD ALWAYS BE 1!
    X_nominal = ekf.X_input(1:ekf.N_states); % This is just our initial condition, not including STM

    [curr_meas, ekf] = EKF.compute_measurement(ekf, S, ENV, X_nominal); % Before updating X

    [X_state_updated, P_cov_updated, dx, ekf] = EKF.ekf_update( ...
        ekf, X_nominal, ekf.P_cov, curr_meas, S, ENV);
    
    EKF.log_step(ekf,dx,P_cov_updated,i);   % LOG before we overwrite ekf

    STM_reset = eye(ekf.N_states);
    ekf.X_input = [X_state_updated; STM_reset(:)]; % Prep for next iteration
    ekf.P_cov   = P_cov_updated;             % Prep for next iteration
    
end