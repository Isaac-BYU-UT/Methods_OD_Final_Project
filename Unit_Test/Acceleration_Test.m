classdef Acceleration_Test < matlab.unittest.TestCase
    methods(Test)

        function test_Accleration_Components(testCase)
            [S,ENV] = Setup.loadSettings('F','Accel',false,false);
        
            time_jd_days = ENV.time_struct_epoch.jd_UTC_days;
            r_ECI_m = S.IC_Sat_Epoch.position_ECI_meters;
            v_ECI_m_s = S.IC_Sat_Epoch.velocity_ECI_meters_per_second;
            C_drag = 1.88;

            % ---- Find Moon and Sun Positions -----
            [r_sun_rel_earth_ECI_meters, ~] = Forces.Vallado_sunPositionLowPrecision(time_jd_days);
            [r_moon_rel_earth_ECI_meters, ~] = Forces.Vallado_moonPositionLowPrecision(time_jd_days);
            sat_is_illuminated = Forces.Vallado_sunLOS(r_ECI_m,r_sun_rel_earth_ECI_meters);
            
            % Compute Rotation Matrix
            date_time_t = datetime(time_jd_days, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
            EOP_params_t = Tools.interpolate_EOP(date_time_t,ENV.EOP_IERS, ENV.EOP_Celestrak);
            [~, ~, R_ECI_from_ECEF, R_ECEF_from_ECI] = Tools.ECEF_ECI_Converter(r_ECI_m, v_ECI_m_s, date_time_t, "ECI_to_ECEF", EOP_params_t);
    
            a_component_ECI_m_s2 = Forces.Compute_Component_Acceleration_ECI_m_s2(...
                                r_ECI_m,v_ECI_m_s,C_drag,...
                                r_sun_rel_earth_ECI_meters(:),...
                                r_moon_rel_earth_ECI_meters(:),...
                                sat_is_illuminated,...
                                R_ECEF_from_ECI);

            a_total_ECI_m_s2 = Forces.Compute_Total_Acceleration_ECI_m_s2(...
                                r_ECI_m,v_ECI_m_s,C_drag,...
                                r_sun_rel_earth_ECI_meters(:),...
                                r_moon_rel_earth_ECI_meters(:),...
                                sat_is_illuminated,...
                                R_ECEF_from_ECI);
    
            % These are all in m/s^2, convert to km/s^2 for comparison
            a_total_ECI_km_s2 = a_total_ECI_m_s2 / Units.KILOMETERS;
            a_component_ECI_km_s2 = a_component_ECI_m_s2 / Units.KILOMETERS;

            disp("a_total_ECI_km_s2: " + mat2str(a_total_ECI_km_s2));
            disp("a_2B: " + mat2str(a_component_ECI_km_s2(1:3)));
            disp("a_Zonals: " + mat2str(a_component_ECI_km_s2(4:6)));
            disp("acc_GRCF_grav:"  + mat2str((a_component_ECI_km_s2(1:3) + a_component_ECI_km_s2(4:6))));
            disp("a_Drag: " + mat2str(a_component_ECI_km_s2(7:9)));
            disp("a_LuniSolar: " + mat2str(a_component_ECI_km_s2(10:12)));
            disp("a_SRP: " + mat2str(a_component_ECI_km_s2(13:15)));


            acc_GRCF_grav_expected = [  -0.0075601602480519;
                                        -0.00175474460199184;
                                        -1.71462729023098e-05];
            acc_drag_expected = [       1.28565567565462e-11;
                                        -5.57684843178588e-11;
                                        -2.15959896971973e-12];
            acc_lunisolar_expected = [     1.03767688038019e-10;
                                            6.44383775293363e-10;
                                            2.98396598163909e-10];
            acc_SRP_expected = [     -2.68799817041286e-11;
                                    -2.9509733168165e-12;
                                    -4.95076071757261e-13];
            a_total_expected = acc_GRCF_grav_expected + acc_drag_expected + acc_lunisolar_expected + acc_SRP_expected;
            disp("Comparison in km/s^2:");
            testCase.verifyEqual(a_component_ECI_km_s2(1:3), acc_GRCF_grav_expected, 'AbsTol', 1e-5);
            testCase.verifyEqual(a_component_ECI_km_s2(7:9), acc_drag_expected, 'RelTol', 1);
            testCase.verifyEqual(a_component_ECI_km_s2(10:12), acc_lunisolar_expected, 'RelTol', 1);
            testCase.verifyEqual(a_component_ECI_km_s2(13:15), acc_SRP_expected, 'RelTol', 1);
            testCase.verifyEqual(a_total_ECI_km_s2, a_total_expected, 'AbsTol', 1e-5);
    
                


        end
    end
end