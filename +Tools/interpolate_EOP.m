function EOP_out = interpolate_EOP(UTC_date_time, EOP_IERS, EOP_Celestrak)

    EOP_Baseline = EOP_IERS;

    % Convert time once
    time_struct = Tools.ComputeTimeSystems(UTC_date_time);
    t_query = time_struct.mjd_UTC_days;

    % Pull pre-packed data (FAST local references)
    t = EOP_Baseline.t;
    Y = EOP_Baseline.Y;

    % ----------------------------
    % FAST index search (linear time if used, but MATLAB optimized)
    % ----------------------------
    i = find(t <= t_query, 1, 'last');

    % Boundary safety (important for EKF robustness)
    if i >= length(t)
        i = length(t) - 1;
    elseif i < 1
        i = 1;
    end

    t0 = t(i);
    t1 = t(i+1);

    % Avoid divide-by-zero (rare but safe)
    if t1 == t0
        alpha = 0;
    else
        alpha = (t_query - t0) / (t1 - t0);
    end

    % ----------------------------
    % LINEAR INTERPOLATION (vectorized row)
    % ----------------------------
    interp_matrix = Y(i,:) + alpha .* (Y(i+1,:) - Y(i,:));

    % ----------------------------
    % Pack output
    % ----------------------------
    EOP_out = struct();

    EOP_out.MJD_days = t_query;

    EOP_out.x_pole_arcsec          = interp_matrix(1);
    EOP_out.y_pole_arcsec          = interp_matrix(2);
    EOP_out.dPsi_milli_arcsec      = interp_matrix(3);
    EOP_out.dEpsilon_milli_arcsec  = interp_matrix(4);
    EOP_out.UT1_UTC_sec            = interp_matrix(5);
    EOP_out.LOD_millisec           = interp_matrix(6);

    % ----------------------------
    % Secondary derived quantities
    % ----------------------------

    EOP_out.delta_AT_sec = interp1( ...
        EOP_Celestrak.MJD_days, ...
        EOP_Celestrak.TAI_minus_UTC_sec, ...
        t_query, ...
        'linear');

    EOP_out.omega_earth_rad_sec = PhysicsConstants.OMEGA_EARTH_RAD_S * ...
        (1 - (EOP_out.LOD_millisec * Units.MILLI_TO_NOM / Units.SEC_IN_SOLAR_DAY));

    EOP_out.small_d_delta_psi_1980_deg = ...
        EOP_out.dPsi_milli_arcsec * Units.MILLI_TO_NOM * Units.ARCSEC_TO_DEG;

    EOP_out.small_d_delta_epsilon_1980_deg = ...
        EOP_out.dEpsilon_milli_arcsec * Units.MILLI_TO_NOM * Units.ARCSEC_TO_DEG;

end