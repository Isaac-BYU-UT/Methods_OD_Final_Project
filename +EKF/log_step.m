function log_step(ekf, dx, P_cov_updated, i)

    if ekf.print_updates && mod(i,ekf.f_updates)==0
        fprintf('Step %d / %d\n', i, ekf.N_obs);
        disp('State correction:'); disp(transpose(dx));
        disp('Covariance Update:'); disp(P_cov_updated);
        % disp('Covariance Difference:'); disp(P_cov_updated - ekf.P_cov)
    end
end