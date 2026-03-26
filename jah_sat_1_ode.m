function dXdt = jah_sat_1_ode(t, X, accel_func,A_func,epoch_jd_UTC_days)

    r = X(1:3);
    v = X(4:6);
    C_drag = X(7);
    STM = reshape(X(8:56),7,7); % 7 states, therefor this will be a 7x7 STM

    % Convert time to Julian Date in days
    time_jd_days = t/Constants.SEC_IN_SOLAR_DAY + epoch_jd_UTC_days;

    % % Compute sun and moon positions at time
    % [r_sun_rel_earth_ECI_km, ~] = planetEphemeris(time_jd_days, 'Earth', 'Sun');
    % [r_moon_rel_earth_ECI_km, ~] = planetEphemeris(time_jd_days, 'Earth', 'Moon');

    r_sun_rel_earth_ECI_km = zeros(1,3);
    r_moon_rel_earth_ECI_km = zeros(1,3);

    A = A_func(r,v,C_drag,r_sun_rel_earth_ECI_km', r_moon_rel_earth_ECI_km'); % Everything must be 3x1

    STM_dot = A*STM; % Propogate STM Matrix
    
    a_total = accel_func(r,v,C_drag,r_sun_rel_earth_ECI_km', r_moon_rel_earth_ECI_km');

    C_drag_dot = 0; % The time derivative of the drag coefficient should be 0.

    dXdt = [v; a_total; C_drag_dot; STM_dot(:)];
end