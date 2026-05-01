function a_total = compute_total_acceleration(t, X, S, ENV)
    % COMPUTE_TOTAL_ACCELERATION Computes total acceleration from current state.
    % This is called by get_A_matrix_FD for finite difference Jacobian calculation.
    %
    % Extracts the acceleration computation from jah_sat_1_ode.m for reuse.
    
    r_ECI_meters = X(1:3);
    v_ECI_meters_s = X(4:6);
    C_drag = SatelliteProperties.C_Drag;
    
    % Convert time to Julian Date in days
    time_jd_days = t / Units.SEC_IN_SOLAR_DAY + ENV.time_struct_epoch.jd_UTC_days;
    
    % ---- Find Moon and Sun Positions -----
    [r_sun_rel_earth_MOD_meters, ~] = Forces.Vallado_sunPositionLowPrecision(time_jd_days);
    [r_moon_rel_earth_MOD_meters, ~] = Forces.Vallado_moonPositionLowPrecision(time_jd_days);
    sat_is_illuminated = Forces.Vallado_sunLOS(r_ECI_meters, r_sun_rel_earth_MOD_meters);
    
    % Compute Rotation Matrix
    date_time_t = datetime(time_jd_days, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
    EOP_params_t = Tools.interpolate_EOP(date_time_t, ENV.EOP_IERS, ENV.EOP_Celestrak);
    [~, ~, ~, R_ECEF_from_ECI, P_Matrix] = Tools.ECEF_ECI_Converter(...
        r_ECI_meters, v_ECI_meters_s, date_time_t, "ECI_to_ECEF", EOP_params_t);
    
    r_sun_rel_earth_ECI_meters = P_Matrix * r_sun_rel_earth_MOD_meters;
    r_moon_rel_earth_ECI_meters = P_Matrix * r_moon_rel_earth_MOD_meters;
    
    % ---- Compute Accelerations ----
    a_NonGravs_m_s2 = Forces.Compute_Total_Acceleration_NonGravs_ECI_m_s2(...
        r_ECI_meters, v_ECI_meters_s, C_drag, ...
        r_sun_rel_earth_ECI_meters(:), ...
        r_moon_rel_earth_ECI_meters(:), ...
        sat_is_illuminated, ...
        R_ECEF_from_ECI);
    
    % EGM96 20x20 spherical harmonic gravity
    L_max = ENV.EGM96_20_x_20_Data.L_max;
    S_coef_matrix = ENV.EGM96_20_x_20_Data.S_coeff_matrix;
    C_coef_matrix = ENV.EGM96_20_x_20_Data.C_coeff_matrix;
    
    mu_m3_s2 = PhysicsConstants.MU_EARTH_KM3_S2 * (Units.KILOMETERS ^ 3);
    re_m = PhysicsConstants.R_EARTH_KM * Units.KILOMETERS;
    
    a_Spherical_m_s2 = Forces.gottliebnorm(mu_m3_s2, re_m, r_ECI_meters, ...
        C_coef_matrix, S_coef_matrix, L_max, L_max, R_ECEF_from_ECI);
    
    % Total acceleration
    a_total = a_Spherical_m_s2 + a_NonGravs_m_s2;
    
end
