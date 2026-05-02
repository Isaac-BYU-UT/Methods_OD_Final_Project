function [y_lt_final_guess_m]  = Light_Time_Correction(initial_range_estimate, ...
                                                            r_sat_original_ECI_m,...
                                                            v_sat_original_ECI_m_s,...
                                                            ekf, S, ENV, curr_meas)

i = ekf.current_index;
last_sat_estimate_ECI_m = r_sat_original_ECI_m;
lt_guess = initial_range_estimate/PhysicsConstants.C_SPEED_OF_LIGHT_M_S;

STM_Dummy = NaN(ekf.N_states,ekf.N_states);
y0_lt = [r_sat_original_ECI_m; v_sat_original_ECI_m_s; STM_Dummy(:)];

sat_delta_m = 999; % initialize large for now
tolerance_m = 1 * Units.MILIMETERS;
num_iterations = 0;

while (any(sat_delta_m > tolerance_m))
    num_iterations = num_iterations + 1;
    t = ekf.t_obs(i); % What's the time step right now
    t_lt = t - lt_guess;
    % JD_lt = ekf.time_struct_epoch.jd_UTC_days + 
    

    if strcmp(ekf.ode_type, 'ode113')

        [~, y_ode_lt] = ode113(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, ekf.debug_on), ...
        [t, t_lt], y0_lt, ekf.options); % ALWAYS PROPOGATE FROM SAME INITIAL CONDITION

    else

        [~, y_ode_lt] = ode45(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, ekf.debug_on), ...
        [t, t_lt], y0_lt, ekf.options); % ALWAYS PROPOGATE FROM SAME INITIAL CONDITION

    end
    
    y_ode_lt_last_state = transpose(y_ode_lt(end,:));
    X_lt = y_ode_lt_last_state(1:ekf.N_states);

    lt_r_sat_ECI_m = X_lt(1:3);
    lt_v_sat_ECI_m_s = X_lt(4:6);
    % Use X_lt to now compute the range , but we need to convert
    % station coordinates to ECI at JD_lt

                                                                    
    date_time_at_lt = S.IC_Sat_Epoch.epoch_date_time_UTC + seconds(t_lt); % Subtract lt_guess!

    EOP_at_lt = Tools.interpolate_EOP(date_time_at_lt, ENV.EOP_IERS, ENV.EOP_Celestrak);

    % Station state in ECI for lt state
    [lt_r_stn_ECI_m, lt_v_stn_ECI_m_s] = Tools.ECEF_ECI_Converter( ...
                                                S.Stations(curr_meas.station_id).position_ECEF_meters, ...
                                                zeros(3,1), date_time_at_lt, "ECEF_to_ECI", EOP_at_lt);

    y_lt_guess_m = Measurements.Compute_Range_Range_Rate(lt_r_sat_ECI_m, lt_v_sat_ECI_m_s,...
                                                lt_r_stn_ECI_m, lt_v_stn_ECI_m_s);

    range_lt_guess_m = y_lt_guess_m(1);
    sat_delta_m = abs(lt_r_sat_ECI_m - last_sat_estimate_ECI_m);
    last_sat_estimate_ECI_m = lt_r_sat_ECI_m;
    lt_guess = range_lt_guess_m/PhysicsConstants.C_SPEED_OF_LIGHT_M_S;

end

if (mod(i, ekf.f_updates) == 0)
    fprintf("Num Iterations: %d ", num_iterations);
end

y_lt_final_guess_m = y_lt_guess_m;


end