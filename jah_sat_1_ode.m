function X_state_dot = jah_sat_1_ode(t, X, S, ENV, debug_on)

    if nargin < 4
        debug_on = false;
    end
    
    r_ECI_meters = X(1:3);
    v_ECI_meters_s = X(4:6);
    C_drag = SatelliteProperties.C_Drag;
    STM = reshape(X(7:42),6,6); % Now only 6 states! % 7 states, therefor this will be a 7x7 STM

    % % Convert time to Julian Date in days
    % time_jd_days = t/Units.SEC_IN_SOLAR_DAY + ENV.time_struct_epoch.jd_UTC_days;

    % % ---- Find Moon and Sun Positions -----
    % [r_sun_rel_earth_MOD_meters, ~] = Forces.Vallado_sunPositionLowPrecision(time_jd_days);
    % [r_moon_rel_earth_MOD_meters, ~] = Forces.Vallado_moonPositionLowPrecision(time_jd_days);
    % sat_is_illuminated = Forces.Vallado_sunLOS(r_ECI_meters,r_sun_rel_earth_MOD_meters);

    %  % Compute Rotation Matrix
    % date_time_t = datetime(time_jd_days, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
    % EOP_params_t = Tools.interpolate_EOP(date_time_t,ENV.EOP_IERS, ENV.EOP_Celestrak);
    % [~, ~, ~, R_ECEF_from_ECI,P_Matrix] = Tools.ECEF_ECI_Converter(r_ECI_meters, v_ECI_meters_s, date_time_t, "ECI_to_ECEF", EOP_params_t);


    % r_sun_rel_earth_ECI_meters = P_Matrix * r_sun_rel_earth_MOD_meters;
    % r_moon_rel_earth_ECI_meters = P_Matrix * r_moon_rel_earth_MOD_meters;


    % Compute and Propogate STM
    % -------------------------
    % Using finite difference A-matrix to avoid 20+ min symbolic compilation
    A_matrix = Forces.get_A_matrix_FD(t, X, S, ENV);
    
    STM_dot = A_matrix*STM; 

    % if t == 0
    %     disp("A_matrix: "); disp((A_matrix - S.ref_data.A_t0_ref)./S.ref_data.A_t0_ref);
    % end

    % --- Compute State Derivatives ---
    % a_total_ECI_m_s2 = Forces.Compute_Total_Acceleration_ECI_m_s2(...
    %                                 r_ECI_m,v_ECI_m_s,C_drag,...
    %                                 r_sun_rel_earth_ECI_meters(:),...
    %                                 r_moon_rel_earth_ECI_meters(:),...
    %                                 sat_is_illuminated,...
    %                                 R_ECEF_from_ECI);

    % a_total_ECI_m_s2 = Forces.Compute_Total_Acceleration_Complex_ECI_m_s2(...
    %                             r_ECI_m,v_ECI_m_s,C_drag,...
    %                             r_sun_rel_earth_ECI_meters(:),...
    %                             r_moon_rel_earth_ECI_meters(:),...
    %                             sat_is_illuminated,...
    %                             R_ECEF_from_ECI,...
    %                             C_coef_matrix, S_coef_matrix, P_coef_matrix);

    % a_2B        = Forces.Gravity_2Body(r_ECI_meters);

    % a_Zonals    = Forces.Gravity_Zonal(r_ECI_meters,...
    %                                     R_ECEF_from_ECI,...
    %                                     true, true, true); % J2, J3, J4 toggles.
    
    % a_Drag      = Forces.Atmospheric_Drag(r_ECI_meters,...
    %                                       v_ECI_meters_s,...
    %                                       C_drag, ...
    %                                       r_sun_rel_earth_ECI_meters);
    
    % a_LuniSolar = Forces.Luni_Solar_Pertubations(r_ECI_meters, ...
    %                                              r_sun_rel_earth_ECI_meters,...
    %                                              r_moon_rel_earth_ECI_meters);
    
    % a_SRP       = Forces.Solar_Radiation_Pressure(r_ECI_meters,...
    %                                               r_sun_rel_earth_ECI_meters,...
    %                                               sat_is_illuminated);
    
    % a_Spherical = Forces.Gravity_Spherical(r_ECI_meters,...
    %                                        R_ECEF_from_ECI,...
    %                                        C_coef_matrix,...
    %                                        S_coef_matrix,...
    %                                        P_coef_matrix);
    

    a_total_ECI_m_s2 = Forces.compute_total_acceleration(t, X, S, ENV);

    % C_drag_dot = 0;

    X_state_dot =  [v_ECI_meters_s;...
                    a_total_ECI_m_s2;...
                    % C_drag_dot;...
                    STM_dot(:)];


    % For Debug Purposes Only
    % ------------------------
    % if (debug_on && t == 0.0)
    %     a_component_ECI_m_s2 = Forces.Compute_Component_Acceleration_ECI_m_s2(...
    %                                     r_ECI_m,v_ECI_m_s,C_drag,...
    %                                     r_sun_rel_earth_ECI_meters(:),...
    %                                     r_moon_rel_earth_ECI_meters(:),...
    %                                     sat_is_illuminated,...
    %                                     R_ECEF_from_ECI);
    % 
    %     disp("a_total_ECI_m_s2: " + mat2str(a_total_ECI_m_s2/Units.KILOMETERS));
    %     disp("a_2B: " + mat2str(a_component_ECI_m_s2(1:3)/Units.KILOMETERS));
    %     disp("a_Zonals: " + mat2str(a_component_ECI_m_s2(4:6)/Units.KILOMETERS));
    %     disp("acc_GRCF_grav:"  + mat2str((a_component_ECI_m_s2(1:3) + a_component_ECI_m_s2(4:6))/Units.KILOMETERS))
    %     disp("a_Drag: " + mat2str(a_component_ECI_m_s2(7:9)/Units.KILOMETERS));
    %     disp("a_LuniSolar: " + mat2str(a_component_ECI_m_s2(10:12)/Units.KILOMETERS));
    %     disp("a_SRP: " + mat2str(a_component_ECI_m_s2(13:15)/Units.KILOMETERS));
    % 
    %     % This will pause execution until you press any key in the Command Window
    %     disp('Press any key to continue to the next observation...');
    %     pause;
    % end

end
