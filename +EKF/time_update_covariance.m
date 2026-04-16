function P_bar = time_update_covariance(ekf, STM_propogated)

    Q = zeros(ekf.N_states); % process noise OFF

    P_bar = STM_propogated * ekf.P_cov * transpose(STM_propogated) + Q;
    % ekf.P_cov holds the previous P_cov_i_minus_1
end