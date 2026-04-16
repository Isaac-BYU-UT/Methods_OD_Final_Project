function meas = compute_measurement(ekf, S, ENV, X_nominal, i)

    meas.station_id = S.ref_data.Actual_Measurements.station_id(i);

    meas.y_meas = S.ref_data.Actual_Measurements{i, ...
        {'apparent_range_km','apparent_range_rate_km_s'}}(:) * Units.KILOMETERS; % CONVERT TO METERS

    meas.R = S.Stations(meas.station_id).Covariance;

    % Time
    date_time = ekf.time_struct.jd_UTC_days + seconds(ekf.t_obs(i));

    EOP = Tools.interpolate_EOP(date_time, ENV.EOP_IERS, ENV.EOP_Celestrak);

    % Station state
    [meas.r_stn_ECI_m, meas.v_stn_ECI_m_s] = Tools.ECEF_ECI_Converter( ...
        S.Stations(meas.station_id).position_ECEF_meters, ...
        zeros(3,1), date_time, "ECEF_to_ECI", EOP);

    % Prefit
    meas.r_sat_ECI_m = X_nomial(1:3);
    meas.v_sat_ECI_m_s = X_nominal(4:6);

    meas.y_computed = Measurements.Compute_Range_Range_Rate( ...
        meas.r_sat_ECI_m, meas.v_sat_ECI_m_s, meas.r_stn_ECI_m, meas.v_stn_ECI_m_s);

    meas.residual = meas.y_meas - meas.y_computed; % In meters and meters per second

    ekf.Y_prefit(:,i) = meas.y_computed;
end