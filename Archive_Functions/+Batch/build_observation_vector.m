function [y_batch, curr_meas_all] = build_observation_vector(batch, S, ENV, X_states_all, curr_meas_all)
% Build observation vector for batch least squares
% Includes light-time corrections and computed measurements

    N_obs = batch.N_obs;
    
    % Pre-allocate observation vector (2*N_obs for range and range-rate)
    y_batch = zeros(2 * N_obs, 1);
    y_computed = zeros(2 * N_obs, 1);

    % Build observation vector element by element
    for i = 1:N_obs
    
        batch.current_index = i;

        station_id = curr_meas_all(i).station_id;

        % Get state at this time
        r_sat_ECI_m = X_states_all(1:3, i);
        v_sat_ECI_m_s = X_states_all(4:6, i);

        % Compute measurements without light-time correction first
        y_computed_no_lt = Measurements.Compute_Range_Range_Rate( ...
            r_sat_ECI_m, v_sat_ECI_m_s, ...
            curr_meas_all(i).r_stn_ECI_m, curr_meas_all(i).v_stn_ECI_m_s);

        % Apply light-time correction
        initial_range_estimate = y_computed_no_lt(1);
        [y_computed_lt] = Measurements.Light_Time_Correction( ...
            initial_range_estimate, ...
            r_sat_ECI_m, v_sat_ECI_m_s, ...
            batch, S, ENV, curr_meas_all(i));

        % Apply station bias (Arecibo example)
        if (station_id == 3)
            y_computed_lt(1) = y_computed_lt(1) + S.Arecibo_Range_Bias_m;
        end

        % Store computed measurement
        curr_meas_all(i).y_computed = y_computed_lt;

        % Build observation vector (O-C form: observation minus computed)
        y_batch(2*i-1:2*i, 1) = curr_meas_all(i).y_obs_meters - y_computed_lt;
        y_computed(2*i-1:2*i, 1) = y_computed_lt;

        % Store prefit values in batch structure for plotting
        batch.Y_obs(:, i) = curr_meas_all(i).y_obs_meters;
        batch.Y_computed(:, i) = y_computed_lt;
    end

end
