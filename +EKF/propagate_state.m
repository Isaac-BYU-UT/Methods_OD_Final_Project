function [X_states_propogated, STM_propogated, y_ode_prop, delta_t] = propagate_state(ekf,S, ENV)

    i = ekf.current_index;

    if strcmp(ekf.ode_type, 'ode113')

        [~, y_ode_prop] = ode113(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, ekf.debug_on), ...
        [ekf.t_obs(i-1), ekf.t_obs(i)], ekf.X_input, ekf.options);

    else
        [~, y_ode_prop] = ode45(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, ekf.debug_on), ...
        [ekf.t_obs(i-1), ekf.t_obs(i)], ekf.X_input, ekf.options);

    end

    delta_t = ekf.t_obs(i) - ekf.t_obs(i-1);

    X_full_propagated = transpose(y_ode_prop(end,:)); % Make this a column vector

    X_states_propogated = X_full_propagated(1:ekf.N_states);

    STM_propogated = reshape(X_full_propagated(ekf.N_states+1:end), ...
                  ekf.N_states, ekf.N_states);
end