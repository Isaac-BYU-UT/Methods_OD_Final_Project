function batch = run_batch(S, ENV, ekf)

%% SETTINGS
max_iters = 10;
tol = 1e-6;

batch.X0 = ekf.X_input(1:ekf.N_states);   % Initial guess
batch.P0 = ekf.P_cov;

for iter = 1:max_iters

    fprintf('\nBatch Iteration %d\n', iter);

    % --- Initialize normal equations ---
    A = zeros(ekf.N_states);
    N = zeros(ekf.N_states,1);

    % Reset trajectory
    STM0 = eye(ekf.N_states);
    X_ref = batch.X0;

    t_prev = ekf.t_obs(1);

    for i = 1:ekf.N_obs

        t_i = ekf.t_obs(i);

        % ---- PROPAGATE STATE + STM ----
        X_input = [X_ref; reshape(STM0,[],1)];

        if strcmp(ekf.ode_type,'ode113')
            [~, y] = ode113(@(t,X) jah_sat_1_ode(t,X,S,ENV,false), ...
                           [t_prev, t_i], X_input, ekf.options);
        else
            [~, y] = ode45(@(t,X) jah_sat_1_ode(t,X,S,ENV,false), ...
                           [t_prev, t_i], X_input, ekf.options);
        end

        X_full = y(end,:)';
        X_ref  = X_full(1:ekf.N_states);

        STM = reshape(X_full(ekf.N_states+1:end), ...
                      ekf.N_states, ekf.N_states);

        % ---- MEASUREMENT ----
        [curr_meas, ekf] = EKF.compute_measurement( ...
            ekf, S, ENV, X_ref);

        if ~ekf.station_meas_mask(curr_meas.station_id)
            t_prev = t_i;
            STM0 = STM;
            continue;
        end

        % Residual
        y_i = curr_meas.y - curr_meas.y_hat;

        % Measurement Jacobian (w.r.t current state)
        H_t = curr_meas.H;   % you already compute this in EKF

        % Map to initial state
        H_i = H_t * STM;

        % Weighting
        R_i = curr_meas.R;

        % ---- ACCUMULATE NORMAL EQUATIONS ----
        A = A + H_i' / R_i * H_i;
        N = N + H_i' / R_i * y_i;

        % Prepare next step
        t_prev = t_i;
        STM0 = STM;

    end

    % ---- SOLVE NORMAL EQUATIONS ----
    dx0 = A \ N;

    fprintf('||dx0|| = %.3e\n', norm(dx0));

    % ---- UPDATE INITIAL STATE ----
    batch.X0 = batch.X0 + dx0;

    % ---- CHECK CONVERGENCE ----
    if norm(dx0) < tol
        disp('Batch converged.');
        break;
    end

end

% Covariance
batch.P0 = inv(A);

end