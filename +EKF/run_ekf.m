clear; clc;

%% SETUP
[S, ENV] = Setup.loadSettings('F','1 Hour','Final - 1 Day Subset',true,true);

ekf = initialize_ekf(S, ENV);

%% FIRST OBS
if ekf.t_obs(1) == 0
    ekf = process_first_observation(ekf, S, ENV);
end

%% MAIN LOOP
for i = 2:ekf.N_obs

    ekf.current_index = i;

    % ---- PROPAGATE ----
    [X_nominal, STM] = propagate_state(ekf);

    % ---- TIME UPDATE ----
    P_bar = time_update_covariance(ekf, STM);

    % ---- MEASUREMENT ----
    meas = compute_measurement(ekf, S, ENV, X_nominal, i);

    % ---- UPDATE ----
    [X_nominal, P_cov, dx] = ekf_update( ...
        ekf, X_nominal, P_bar, meas, S, ENV);

    % ---- POSTFIT ----
    ekf.Y_postfit(:,i) = compute_observation( ...
        X_nominal, meas.r_stn, meas.v_stn);

    % ---- LOG ----
    log_step(ekf, dx, P_cov, i);

    % ---- PREP NEXT ----
    STM_reset = reshape(eye(ekf.N_states),[],1);
    ekf.X_input = [X_nominal; STM_reset];
    ekf.P_cov   = P_cov;
end