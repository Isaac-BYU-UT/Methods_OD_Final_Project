function batch = log_batch(batch, dX_batch, P_cov_batch, varargin)
% Logging function for batch least squares OD
% Following similar conventions as EKF.log_step

    if batch.print_updates
        fprintf('=== Batch Least Squares Solution ===\n');
        fprintf('Number of observations: %d\n', batch.N_obs);
        fprintf('Number of states: %d\n', batch.N_states);
        fprintf('Degrees of freedom: %d\n', 2*batch.N_obs - batch.N_states);
        fprintf('\nState Correction (dX):\n');
        disp(transpose(dX_batch));
        fprintf('\nState Correction Covariance (diagonal):\n');
        disp(diag(P_cov_batch)');
        fprintf('=========================================\n\n');
    end

    % Optional: Compute and display chi-squared statistics if weight matrix provided
    if ~isempty(varargin) && ~isempty(varargin{1})
        W = varargin{1};
        y = varargin{2};
        residuals = y - (varargin{3} * dX_batch);  % H * dX
        
        chi_sq = transpose(residuals) * W * residuals;
        dof = 2*batch.N_obs - batch.N_states;
        
        if batch.print_updates
            fprintf('Chi-squared statistic: %e\n', chi_sq);
            fprintf('Normalized chi-squared (chi-sq/DOF): %e\n', chi_sq/dof);
        end
    end

end
