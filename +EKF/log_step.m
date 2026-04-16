function log_step(ekf, dx, P_cov, i)

    if evalin('base','print_updates') || mod(i,20)==0
        fprintf('Step %d / %d\n', i, ekf.N_obs);
        disp('State correction:'); disp(transpose(dx));
        disp('Covariance:'); disp(P_cov);
    end
end