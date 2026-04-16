classdef ECEF_ECI_Test < matlab.unittest.TestCase
    methods(Test)
        % function test_ECI_to_ECEF_Vallado(testCase)
        %     %% Postiion and Velocity Vectors
        %     % Radius X Y Z ECEF (ITRF) [in kilometers]-1033.4793830 7901.2952754 6380.3565958
        %     % Velocity X Y Z ECEF (ITRF) [in kilometers]-3.225636520-2.872451450 5.531924446
        % 
        %     r_ECEF = [-1033.4793830, 7901.2952754, 6380.3565958]';
        %     v_ECEF = [-3.225636520, -2.872451450, 5.531924446]';
        % 
        %     % Radius X Y Z ECI (GCRF) [in kilometers] 5102.508958 6123.011401 6378.136928
        %     % Velocity X Y Z ECI (GCRF) [in kilometers]-4.74322016 0.79053650 5.53375528
        % 
        %     r_ECI = [5102.508958 6123.011401 6378.136928]';
        %     v_ECI = [-4.74322016 0.79053650 5.53375528]';
        % 
        %     % Gregorian Date (UTC) Year: 2004 Month: April Day: 6 Hour: 7 Minute: 51 Seconds: 28.386009
        %     UTC.year = 2004;
        %     UTC.month = 4;
        %     UTC.day = 6;
        %     UTC.hour = 7;
        %     UTC.minute = 51;
        %     UTC.seconds = 28.386009;
        % 
        %     UTC_date_time = datetime(UTC.year, UTC.month, UTC.day, ...
        %                        UTC.hour, UTC.minute, UTC.seconds, ...
        %                        'TimeZone','UTC');
        % 
        %     EOP_params = Tools.interpolate_EOP(UTC_date_time,"IERS");
        % 
        % 
        %     [r_ECI_actual,v_ECI_actual,~,~] = Tools.ECEF_ECI_Converter(r_ECEF,v_ECEF,UTC_date_time,'ECEF_to_ECI',EOP_params);
        % 
        % 
        %     r_ECI_expected = dcmeci2ecef('IAU-76/FK5', UTC_date_time, EOP_params.delta_AT_sec) \ r_ECEF;
        % 
        %     diff_r_km = abs(r_ECI_actual - r_ECI_expected)
        %     r_ECI_expected - r_ECI
        %     r_ECI_actual - r_ECI
        % 
        %     testCase.verifyLessThanOrEqual(diff_r_km,1e-3);
        % 
        % 
        %     [r_ECEF_actual,v_ECEF_actual,~,~] = Tools.ECEF_ECI_Converter(r_ECI,v_ECI,UTC_date_time,'ECI_to_ECEF',EOP_params);
        %     r_ECEF_expected = dcmeci2ecef('IAU-76/FK5', UTC_date_time, EOP_params.delta_AT_sec) * r_ECI;
        % 
        %     diff_r_km = abs(r_ECEF_actual - r_ECEF_expected)
        %     r_ECEF_expected - r_ECEF
        %     r_ECEF_actual - r_ECEF
        % 
        %     testCase.verifyLessThanOrEqual(diff_r_km,1e-3);
        % end

        function test_ECI_to_ECEF_Project(testCase)
            [S,ENV] = Setup.loadSettings('F','1 Hour','Final - Full 3 Days');


            [r_ECEF_actual,v_ECEF_actual,R_ECI_from_ECEF_actual, R_ECEF_from_ECI_actual] = Tools.ECEF_ECI_Converter(S.IC_Sat_Epoch.position_ECI_meters,...
                                                                        S.IC_Sat_Epoch.velocity_ECI_meters_per_second, ...
                                                                        S.IC_Sat_Epoch.epoch_date_time_UTC,...
                                                                        'ECI_to_ECEF',...
                                                                        ENV.EOP_t0);
          
            R_ECI_to_ECEF = dcmeci2ecef('IAU-76/FK5',...
                                            S.IC_Sat_Epoch.epoch_date_time_UTC,...
                                            ENV.EOP_t0.delta_AT_sec, ...
                                            ENV.EOP_t0.UT1_UTC_sec,...
                                            [ENV.EOP_t0.x_pole_arcsec, ENV.EOP_t0.y_pole_arcsec] * Units.ARCSEC_TO_RAD,...
                                            'dNutation',[ENV.EOP_t0.dPsi_milli_arcsec,ENV.EOP_t0.dEpsilon_milli_arcsec] * Units.MILLI_TO_NOM * Units.ARCSEC_TO_RAD...
                                            );
            % 
            % xp_rad = ENV.EOP_t0.x_pole_arcsec * Units.ARCSEC_TO_RAD;
            % yp_rad = ENV.EOP_t0.y_pole_arcsec * Units.ARCSEC_TO_RAD;
            % 
            % % Approximate Version
            % M = [1 0  -xp_rad;
            %      0 1    yp_rad;
            %      xp_rad -yp_rad 1];

            r_ECEF_expected = R_ECI_to_ECEF * S.IC_Sat_Epoch.position_ECI_meters;
            % r_ECEF_expected = eci2ecef(S.IC_Sat_Epoch.epoch_date_time_UTC,...
            %                             S.IC_Sat_Epoch.position_ECI_meters,...
            %                             S.IC_Sat_Epoch.velocity_ECI_meters_per_second,...
            %                             [0;0;0], ...
            %                             'dAT',ENV.EOP_t0.delta_AT_sec,...
            %                             'dUT1', ENV.EOP_t0.UT1_UTC_sec,...
            %                             'pm',[ENV.EOP_t0.x_pole_arcsec, ENV.EOP_t0.y_pole_arcsec] * Units.ARCSEC_TO_DEG ...
            %                             );

            % v_ECEF_expected = 0;

            % Omega_Earth_Vector_Rad_Sec = [0; 0; EOP_params.omega_earth_rad_sec];
            % TODO: Come back to this later!

            diff_r_meters = abs(r_ECEF_actual - r_ECEF_expected);
            disp(diff_r_meters);
            disp('rot matrices')
            R_ECEF_from_ECI_actual - R_ECI_to_ECEF
            
            % diff_v_m_per_s = abs(v_ECEF_actual - v_ECEF_expected)

            testCase.verifyLessThanOrEqual(diff_r_meters,2);
            % testCase.verifyLessThanOrEqual(diff_v_m_per_s, 1e-1);
        end
    end
end