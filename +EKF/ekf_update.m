function [X_updated, P_cov, dx] = ekf_update( ...
    ekf, X_nominal, P_bar, meas, S, ENV)

    update = ekf.station_meas_mask(meas.station_id); % Looks inside this mask to see if the station is on.

    if update

        H_tilde = Measurements.Compute_H_matrix( ...
                meas.r_sat_ECI_m, meas.v_sat_ECI_m_s, meas.r_stn_ECI_m, meas.v_stn_ECI_m_s);

        K = P_bar * transpose(H_tilde) / (H_tilde * P_bar * transpose(H_tilde) + meas.R);

        dx = K * meas.residual;

        X_updated = X_nominal + dx;

        I = eye(ekf.N_states);

        P_cov = (I - K*H_tilde)*P_bar*(I - K*H_tilde)' + K*meas.R*K';
        P_cov = 0.5 * (P_cov + P_cov');

    else
        X_updated = X_nominal;
        P_cov = P_bar;
        dx = zeros(ekf.N_states,1);
    end
end