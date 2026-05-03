function [X_states_all, STM_all, y_ode_all] = propagate_states_batch(batch, S, ENV)
% Propagate state and State Transition Matrices for all observation times
% Returns matrices with states and STMs for each observation epoch

    N_states = batch.N_states;
    N_obs = batch.N_obs;
    t_obs = batch.t_obs;

    % Pre-allocate storage for all states and STMs
    X_states_all = zeros(N_states, N_obs);
    STM_all = zeros(N_states, N_states, N_obs);
    y_ode_all = [];

    % Initial state preparation (with identity STM)
    STM0 = eye(N_states);
    X0_augmented = [batch.X_input; STM0(:)];

    % Propagate from initial epoch to all observation times
    if strcmp(batch.ode_type, 'ode113')
        [~, y_ode_temp] = ode113(@(t, X) jah_sat_1_ode( ...
            t, X, S, ENV, batch.debug_on), ...
            t_obs, X0_augmented, batch.options);
    else
        [~, y_ode_temp] = ode45(@(t, X) jah_sat_1_ode( ...
            t, X, S, ENV, batch.debug_on), ...
            t_obs, X0_augmented, batch.options);
    end

    % Store full ODE output for later use
    y_ode_all = y_ode_temp;

    % Extract states and STMs at each observation time
    for i = 1:N_obs
        X_full = transpose(y_ode_temp(i, :)); % Column vector

        % Extract state
        X_states_all(:, i) = X_full(1:N_states);

        % Extract STM
        STM_all(:, :, i) = reshape(X_full(N_states+1:end), N_states, N_states);
    end

end
