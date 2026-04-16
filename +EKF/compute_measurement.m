function [curr_meas,ekf] = compute_measurement(ekf, S, ENV, X_nominal)
% meas struct with fields:
%       station_id
%       y_obs_meters
%       R
%       r_stn_ECI_m
%       v_stn_ECI_m_s
%       y_computed_meters
%       residual
% -----------------------
    i = ekf.current_index;
    curr_meas.station_id = S.ref_data.Actual_Measurements.station_id(i);

    curr_meas.y_obs_meters = [S.ref_data.Actual_Measurements.apparent_range_km(i); ...
                   S.ref_data.Actual_Measurements.apparent_range_rate_km_s(i)] ...
                   * Units.KILOMETERS; % SUPER IMPORTANT!!!

    curr_meas.R = S.Stations(curr_meas.station_id).Covariance;

    % Time
    date_time = S.IC_Sat_Epoch.epoch_date_time_UTC + seconds(ekf.t_obs(i)); % Yep, we want time at t_obs(i)!

    EOP = Tools.interpolate_EOP(date_time, ENV.EOP_IERS, ENV.EOP_Celestrak);

    % Station state
    [curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s] = ...
                                             Tools.ECEF_ECI_Converter( ...
                                                S.Stations(curr_meas.station_id).position_ECEF_meters, ...
                                                zeros(3,1), date_time, "ECEF_to_ECI", EOP);

    % Prefit
    r_sat_nominal_ECI_m = X_nominal(1:3);
    v_sat_nominal_ECI_m_s = X_nominal(4:6);

    curr_meas.y_computed_meters = Measurements.Compute_Range_Range_Rate( ...
                                                                    r_sat_nominal_ECI_m, v_sat_nominal_ECI_m_s,... % Before update
                                                                    curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s); % Comptued ECI

    ekf.Y_prefit(:,i) = curr_meas.y_computed_meters;
    curr_meas.residual = curr_meas.y_obs_meters - curr_meas.y_computed_meters; % In meters and meters per second
end