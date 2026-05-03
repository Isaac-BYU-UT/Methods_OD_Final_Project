function [batch, curr_meas_all] = compute_residuals_batch(batch, S, ENV, X_states_all, X_states_updated, dX_batch, curr_meas_all)
% Compute prefit and postfit residuals for all measurements
% Prefit: residuals with propagated (nominal) state
% Postfit: residuals with updated state (after correction)

    N_obs = batch.N_obs;

    % Prefit residuals already computed in build_observation_vector
    % Now compute postfit residuals with updated state

    % Apply state correction to all states
    X_states_postfit = X_states_all + dX_batch; % Correction is broadcast to all measurements

    % Compute postfit measurements
    for i = 1:N_obs
        
        batch.current_index = i;

        station_id = curr_meas_all(i).station_id;

        % Get updated state at this time
        r_sat_updated_ECI_m = X_states_postfit(1:3, i);
        v_sat_updated_ECI_m_s = X_states_postfit(4:6, i);

        % Compute postfit measurements without light-time correction
        y_postfit_no_lt = Measurements.Compute_Range_Range_Rate( ...
            r_sat_updated_ECI_m, v_sat_updated_ECI_m_s, ...
            curr_meas_all(i).r_stn_ECI_m, curr_meas_all(i).v_stn_ECI_m_s);

        % Apply light-time correction
        initial_range_estimate = y_postfit_no_lt(1);
        [y_postfit_lt] = Measurements.Light_Time_Correction( ...
            initial_range_estimate, ...
            r_sat_updated_ECI_m, v_sat_updated_ECI_m_s, ...
            batch, S, ENV, curr_meas_all(i));

        % Apply station bias
        if (station_id == 3)
            y_postfit_lt(1) = y_postfit_lt(1) + S.Arecibo_Range_Bias_m;
        end

        % Store postfit measurement
        curr_meas_all(i).y_computed_postfit = y_postfit_lt;

        % Compute postfit residual
        residual_postfit = curr_meas_all(i).y_obs_meters - y_postfit_lt;
        batch.Y_residuals(:, i) = residual_postfit;

        % Store postfit values
        batch.Y_postfit(:, i) = y_postfit_lt;

        % Store measurement covariance for plotting
        batch.P_zz(:, :, i) = curr_meas_all(i).R;
    end

end
