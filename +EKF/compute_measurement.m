function [curr_meas,ekf] = compute_measurement(ekf, S, ENV, X_states_propogated)
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
    date_time_at_meas = S.IC_Sat_Epoch.epoch_date_time_UTC + seconds(ekf.t_obs(i)); % Yep, we want time at t_obs(i)!

    EOP_at_meas = Tools.interpolate_EOP(date_time_at_meas, ENV.EOP_IERS, ENV.EOP_Celestrak);

    % Station state
    [curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s] = ...
                                             Tools.ECEF_ECI_Converter( ...
                                                S.Stations(curr_meas.station_id).position_ECEF_meters, ...
                                                zeros(3,1), date_time_at_meas, "ECEF_to_ECI", EOP_at_meas);

    % Prefit
    r_sat_propogated_ECI_m = X_states_propogated(1:3);
    v_sat_propogated_ECI_m_s = X_states_propogated(4:6);


    % THIS IS NOT YET LIGHT TIME CORRECTED!!!
    curr_meas.y_computed_propogated_no_lt_meters = Measurements.Compute_Range_Range_Rate( ...
                                                                    r_sat_propogated_ECI_m, v_sat_propogated_ECI_m_s,... % Before update
                                                                    curr_meas.r_stn_ECI_m, curr_meas.v_stn_ECI_m_s); % Comptued ECI

end