function P_bar = time_update_covariance(ekf, STM_propogated, delta_t, S)

    dt = delta_t;
    s2x = S.StateCovariances.sigma_x_accel_meters_per_sec2 ^ 2;
    s2y = S.StateCovariances.sigma_y_accel_meters_per_sec2 ^ 2;
    s2z = S.StateCovariances.sigma_z_accel_meters_per_sec2 ^ 2;

    Process_Noise = [ (dt^4/4)*s2x,       0,             0,         (dt^3/2)*s2x,       0,             0;
                           0,        (dt^4/4)*s2y,       0,              0,        (dt^3/2)*s2y,       0;
                           0,             0,        (dt^4/4)*s2z,        0,             0,        (dt^3/2)*s2z;
                      (dt^3/2)*s2x,       0,             0,          (dt^2)*s2x,        0,             0;
                           0,        (dt^3/2)*s2y,       0,              0,         (dt^2)*s2y,       0;
                           0,             0,        (dt^3/2)*s2z,        0,             0,         (dt^2)*s2z ];

    P_bar = STM_propogated * ekf.P_cov * transpose(STM_propogated) + Process_Noise;
    % ekf.P_cov holds the previous P_cov_i_minus_1

    % We need to ensure that P_bar represents higher uncertainty than ekf.P_cov!
    % Uncertainty should grow during the time update step, so we can check if P_bar is greater than P_cov.
    if ((ekf.debug_on) && any(diag(P_bar) < diag(ekf.P_cov)))
        warning('P_bar has smaller diagonal elements than P_cov, which may indicate a decrease in uncertainty. Please check the STM and process noise.');
    end
end