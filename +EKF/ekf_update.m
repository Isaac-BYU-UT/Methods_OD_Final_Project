function [X_updated, P_cov_updated, dx] = ekf_update( ...
    ekf, X_states_propogated, P_bar, meas, S, ENV)

% NO CHANGES TO ekf

    update_this_step = ekf.station_meas_mask(meas.station_id); % Looks inside this mask to see if the station is on.

    if update_this_step

        r_sat_nominal_ECI_m = X_states_propogated(1:3);
        v_sat_nominal_ECI_m_s = X_states_propogated(4:6);

        H_tilde = Measurements.Compute_H_matrix( ...
                                                r_sat_nominal_ECI_m, v_sat_nominal_ECI_m_s,... % Before update
                                                meas.r_stn_ECI_m, meas.v_stn_ECI_m_s); % Pull in the station ECI coordinates we computed in compute_measurement

        K = P_bar * transpose(H_tilde) / (H_tilde * P_bar * transpose(H_tilde) + meas.R); % Covariance for that specific station

        dx = K * meas.residual;

        X_updated = X_states_propogated + dx;

        I = eye(ekf.N_states);

        P_cov_updated = (I - K*H_tilde)*P_bar*transpose(I - K*H_tilde) + K*meas.R*transpose(K);
        P_cov_updated = 0.5 * (P_cov_updated + transpose(P_cov_updated)); % Ensure symmetry

    else
        dx = zeros(ekf.N_states,1);
        X_updated = X_states_propogated; % No correction applied
        P_cov_updated = P_bar; % No updated covariance
    end
end