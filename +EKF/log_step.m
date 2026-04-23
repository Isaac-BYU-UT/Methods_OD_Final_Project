function ekf =  log_step(ekf, dx, P_cov_updated, i, varargin)
    % Inputs:
    %   ekf - EKF structure
    %   dx - state correction
    %   P_cov_updated - updated covariance
    %   i - observation index
    %   varargin{1} - optional P_bar (covariance after propagation, before update)

    if ekf.print_updates && mod(i,ekf.f_updates)==0
        fprintf('Step %d / %d\n', i, ekf.N_obs);
        disp('State correction:'); disp(transpose(dx));
        disp('Covariance Update:'); disp(P_cov_updated);
        disp('Covariance Det:'); disp(det(P_cov_updated));
    end
    
    % Record trace before propagation (from previous iteration)
    ekf.trace_pre_propagation(i) = trace(ekf.P_cov);
    
    % Record trace after update
    ekf.trace_post_update(i) = trace(P_cov_updated);
    
    % If P_bar is provided, record its trace
    if ~isempty(varargin) && ~isempty(varargin{1})
        ekf.trace_post_propagation(i) = trace(varargin{1});
    end
end