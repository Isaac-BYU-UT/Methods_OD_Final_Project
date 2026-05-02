function [X_states_updated, P_cov_updated, dx, ekf] = ekf_update(ekf, X_states_propogated, P_bar, curr_meas, S, ENV)

    i = ekf.current_index;
    update_this_step = ekf.station_meas_mask(curr_meas.station_id); % Looks inside this mask to see if the station is on.

    if update_this_step

        % --- Prefit (Covariance) ---
        r_sat_propogated_ECI_m = X_states_propogated(1:3);
        v_sat_propogated_ECI_m_s = X_states_propogated(4:6);

        H_tilde_prefit = Measurements.Compute_H_matrix( r_sat_propogated_ECI_m, v_sat_propogated_ECI_m_s,... % Before update
                                                            curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s); % Pull in the station ECI coordinates we computed in compute_measurement
                    
        
        P_zz_prefit = H_tilde_prefit * P_bar * transpose(H_tilde_prefit) + curr_meas.R;
        ekf.P_zz_prefit(:,:,i) = P_zz_prefit;

        % --- Prefit w/ Light Time Correction (Measurement) ---
        
        initial_range_estimate_prefit = curr_meas.y_computed_propogated_no_lt_meters(1); % Range

        [y_lt_prefit_final_guess_m] = Measurements.Light_Time_Correction(...
                                                    initial_range_estimate_prefit, ...
                                                    r_sat_propogated_ECI_m,...
                                                    v_sat_propogated_ECI_m_s,...
                                                    ekf, S, ENV, curr_meas);


        curr_meas.y_lt_computed_propogated_meters = y_lt_prefit_final_guess_m; % If we set this equal to the last viable, we should be good!
        ekf.Y_prefit(:,i) = curr_meas.y_lt_computed_propogated_meters;

        curr_meas.residual = curr_meas.y_obs_meters - curr_meas.y_lt_computed_propogated_meters;

        % -- Kalman Gain ---
        K = P_bar * transpose(H_tilde_prefit) / (P_zz_prefit); % Covariance for that specific station

        % --- Covariance Updated ----
        I = eye(ekf.N_states);
        % P_cov_updated = (I - K*H_tilde)*P_bar; % Non-Joseph form
        P_cov_updated = (I - K*H_tilde_prefit)*P_bar*transpose(I - K*H_tilde_prefit) + K*curr_meas.R*transpose(K); % Joseph form
        P_cov_updated = 0.5 * (P_cov_updated + transpose(P_cov_updated)); % Ensure symmetry
        

        % --- State Update ---
        dx = K * curr_meas.residual;

        X_states_updated = X_states_propogated + dx;

        r_sat_filtered_ECI_m = X_states_updated(1:3);
        v_sat_filtered_ECI_m_s = X_states_updated(4:6);

        % ---- POSTFIT (Covariance) ----

        H_tilde_postfit = Measurements.Compute_H_matrix( ...
                                                r_sat_filtered_ECI_m, v_sat_filtered_ECI_m_s,...
                                                curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s);
        
        P_zz_postfit = H_tilde_postfit * P_cov_updated * transpose(H_tilde_postfit) + curr_meas.R;
        ekf.P_zz_postfit(:,:,i) = P_zz_postfit; % Store this

        % --- POSTFIT (Measurement) ---
        
        curr_meas.y_computed_filtered_no_lt_meters = Measurements.Compute_Range_Range_Rate(...
                                                            r_sat_filtered_ECI_m, v_sat_filtered_ECI_m_s,...
                                                            curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s);
        
        initial_range_estimate_postfit = curr_meas.y_computed_filtered_no_lt_meters(1); % Range

        [y_lt_postfit_final_guess_m] = Measurements.Light_Time_Correction(...
                                            initial_range_estimate_postfit, ...
                                            r_sat_filtered_ECI_m,...
                                            v_sat_filtered_ECI_m_s,...
                                            ekf, S, ENV, curr_meas);

        curr_meas.y_lt_computed_filtered_meters = y_lt_postfit_final_guess_m;
        ekf.Y_postfit(:,ekf.current_index) = curr_meas.y_lt_computed_filtered_meters;

    else
        dx = zeros(ekf.N_states,1);
        X_states_updated = X_states_propogated; % No correction applied
        P_cov_updated = P_bar; % No updated covariance
    end
end