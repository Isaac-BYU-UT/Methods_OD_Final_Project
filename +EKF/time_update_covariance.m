function P_bar = time_update_covariance(ekf, STM_propogated, delta_t, S)

    dt = delta_t;
    % TODO: Is this RIC process noise?
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

end