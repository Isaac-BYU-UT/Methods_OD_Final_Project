function X_state_dot = jah_sat_1_ode(t, X, epoch_jd_UTC_days,debug_on)

    if nargin < 4
        debug_on = false;
    end
    
    r_ECI_m = X(1:3);
    v_ECI_m_s = X(4:6);
    C_drag = X(7);
    STM = reshape(X(8:56),7,7); % 7 states, therefor this will be a 7x7 STM

    % Convert time to Julian Date in days
    time_jd_days = t/Constants.SEC_IN_SOLAR_DAY + epoch_jd_UTC_days;

    % ---- Find Moon and Sun Positions -----
    [r_sun_rel_earth_ECI_meters, ~] = Forces.Vallado_sunPositionLowPrecision(time_jd_days);
    [r_moon_rel_earth_ECI_meters, ~] = Forces.Vallado_moonPositionLowPrecision(time_jd_days);
    sat_is_illuminated = Forces.Vallado_sunLOS(r_ECI_m,r_sun_rel_earth_ECI_meters);

    % Compute Rotation Matrix
    date_time_t = datetime(time_jd_days, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
    EOP_params_t = Tools.interpolate_EOP(date_time_t,"IERS");
    [~, ~, R_ECI_from_ECEF, R_ECEF_from_ECI] = Tools.ECEF_ECI_Converter(r_ECI_m, v_ECI_m_s, date_time_t, "ECI_to_ECEF", EOP_params_t);

    % Compute and Propogate STM
    % -------------------------
    A_matrix = Forces.get_A_matrix(...
                                    r_ECI_m,v_ECI_m_s,C_drag,...
                                    r_sun_rel_earth_ECI_meters(:),...
                                    r_moon_rel_earth_ECI_meters(:),...
                                    sat_is_illuminated,...
                                    R_ECEF_from_ECI);
    
    STM_dot = A_matrix*STM; 


    % --- Compute State Derivatives ---
    a_total_ECI_m_s2 = Forces.Compute_Total_Acceleration_ECI_m_s2(...
                                    r_ECI_m,v_ECI_m_s,C_drag,...
                                    r_sun_rel_earth_ECI_meters(:),...
                                    r_moon_rel_earth_ECI_meters(:),...
                                    sat_is_illuminated,...
                                    R_ECEF_from_ECI);

    C_drag_dot = 0;

    X_state_dot =  [v_ECI_m_s;...
                    a_total_ECI_m_s2;...
                    C_drag_dot;...
                    STM_dot(:)];


    % For Debug Purposes Only
    % ------------------------
    if debug_on
        a_component_ECI_m_s2 = Forces.Compute_Component_Acceleration_ECI_m_s2(...
                                        r_ECI_m,v_ECI_m_s,C_drag,...
                                        r_sun_rel_earth_ECI_meters(:),...
                                        r_moon_rel_earth_ECI_meters(:),...
                                        sat_is_illuminated,...
                                        R_ECEF_from_ECI);

        disp("a_total_ECI_m_s2: " + mat2str(a_total_ECI_m_s2/Units.KILOMETERS));
        disp("a_2B: " + mat2str(a_component_ECI_m_s2(1:3)/Units.KILOMETERS));
        disp("a_Zonals: " + mat2str(a_component_ECI_m_s2(4:6)/Units.KILOMETERS));
        disp("acc_GRCF_grav:"  + mat2str((a_component_ECI_m_s2(1:3) + a_component_ECI_m_s2(4:6))/Units.KILOMETERS))
        disp("a_Drag: " + mat2str(a_component_ECI_m_s2(7:9)/Units.KILOMETERS));
        disp("a_LuniSolar: " + mat2str(a_component_ECI_m_s2(10:12)/Units.KILOMETERS));
        disp("a_SRP: " + mat2str(a_component_ECI_m_s2(13:15)/Units.KILOMETERS));
    end

end