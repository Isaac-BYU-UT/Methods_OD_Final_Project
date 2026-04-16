function P_bar = time_update_covariance(ekf, STM)

    Q = zeros(ekf.N_states); % process noise OFF

    P_bar = STM * ekf.P_cov * transpose(STM) + Q;
end