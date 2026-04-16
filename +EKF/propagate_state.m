function [X_nominal, STM] = propagate_state(ekf)

    i = ekf.current_index;

    [~, y] = ode45(@(t,X) jah_sat_1_ode( ...
        t, X, ekf.time_struct.jd_UTC_days, false), ...
        [ekf.t_obs(i-1), ekf.t_obs(i)], ekf.X_input, ekf.options);

    X_full = y(end,:)';

    X_nominal = X_full(1:ekf.N_states);

    STM = reshape(X_full(ekf.N_states+1:end), ...
                  ekf.N_states, ekf.N_states);
end