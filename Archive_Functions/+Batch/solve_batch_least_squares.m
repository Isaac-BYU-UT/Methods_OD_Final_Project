function [dX_batch, P_cov_batch, batch] = solve_batch_least_squares(batch, S, H_batch, y_batch, curr_meas_all)
% Solve batch least squares using normal equations
% dX_batch = (H^T * W * H)^-1 * H^T * W * y
% where W is the weight matrix (inverse of measurement covariance)

    N_obs = batch.N_obs;
    N_states = batch.N_states;

    % Build weight matrix (inverse of measurement covariance)
    % For uncorrelated measurements, this is diagonal
    W_diag = zeros(2 * N_obs, 1);

    for i = 1:N_obs
        R_i = curr_meas_all(i).R;
        R_i_inv = inv(R_i);
        
        W_diag(2*i-1:2*i) = diag(R_i_inv);
    end

    W = diag(W_diag);

    % Normal equations for batch least squares
    % Build H^T * W
    HTW = transpose(H_batch) * W;

    % Build normal matrix (H^T * W * H)
    N_normal = HTW * H_batch;

    % Build weighted observation vector (H^T * W * y)
    HTWy = HTW * y_batch;

    % Check condition number of normal matrix
    cond_num = cond(N_normal);
    if batch.print_updates
        fprintf('Condition number of normal matrix: %e\n', cond_num);
    end

    % Solve for state correction: dX = (H^T * W * H)^-1 * H^T * W * y
    dX_batch = N_normal \ HTWy;

    % Compute covariance of the solution
    % P_cov = (H^T * W * H)^-1
    P_cov_batch = inv(N_normal);

    % Store matrices for diagnostics
    batch.H_full = H_batch;
    batch.W_full = W;
    batch.N_normal = N_normal;

    if batch.print_updates
        fprintf('State correction (dX):\n');
        disp(transpose(dX_batch));
        fprintf('Solution covariance diagonal:\n');
        disp(diag(P_cov_batch)');
    end

end
