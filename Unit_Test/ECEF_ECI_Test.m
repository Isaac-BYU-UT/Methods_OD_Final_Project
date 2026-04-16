classdef ECEF_ECI_Test < matlab.unittest.TestCase
    methods(Test)
        function test_ECI_to_ECEF_Vallado(testCase)
            [S,ENV] = Setup.loadSettings('F','ECI_ECEF_test',false,false);
            %% Postiion and Velocity Vectors
            % Radius X Y Z ECEF (ITRF) [in kilometers]-1033.4793830 7901.2952754 6380.3565958
            % Velocity X Y Z ECEF (ITRF) [in kilometers]-3.225636520-2.872451450 5.531924446

            r_ECEF_ref = [-1033.4793830, 7901.2952754, 6380.3565958]';
            v_ECEF_ref = [-3.225636520, -2.872451450, 5.531924446]';

            % Radius X Y Z ECI (GCRF) [in kilometers] 5102.508958 6123.011401 6378.136928
            % Velocity X Y Z ECI (GCRF) [in kilometers]-4.74322016 0.79053650 5.53375528

            r_ECI_ref = [5102.508958 6123.011401 6378.136928]';
            v_ECI_ref = [-4.74322016 0.79053650 5.53375528]';

            % Gregorian Date (UTC) Year: 2004 Month: April Day: 6 Hour: 7 Minute: 51 Seconds: 28.386009
            UTC.year = 2004;
            UTC.month = 4;
            UTC.day = 6;
            UTC.hour = 7;
            UTC.minute = 51;
            UTC.seconds = 28.386009;

            UTC_date_time = datetime(UTC.year, UTC.month, UTC.day, ...
                               UTC.hour, UTC.minute, UTC.seconds, ...
                               'TimeZone','UTC');

            EOP_params = Tools.interpolate_EOP(UTC_date_time,ENV.EOP_IERS,ENV.EOP_Celestrak);


            [r_ECI_actual,v_ECI_actual,~,~] = Tools.ECEF_ECI_Converter(r_ECEF_ref,v_ECEF_ref,UTC_date_time,'ECEF_to_ECI',EOP_params);


            r_ECI_builtin = dcmeci2ecef('IAU-76/FK5', UTC_date_time, EOP_params.delta_AT_sec) \ r_ECEF_ref;



            [r_ECEF_actual,v_ECEF_actual,~,~] = Tools.ECEF_ECI_Converter(r_ECI_ref,v_ECI_ref,UTC_date_time,'ECI_to_ECEF',EOP_params);
            r_ECEF_expected = dcmeci2ecef('IAU-76/FK5', UTC_date_time, EOP_params.delta_AT_sec) * r_ECI_ref;

            disp('ECI Position Comparison (Actual vs Built-in):');
            disp(r_ECI_actual);
            disp(r_ECI_builtin);
            disp('ECI Position Comparison (Built-in vs Reference):');
            disp(r_ECI_builtin);
            disp(r_ECI_ref);
            disp('ECI Velocity Comparison (Actual vs Reference):');
            disp(v_ECI_actual);
            disp(v_ECI_ref);
            disp('ECEF Position Comparison (Actual vs Expected):');
            disp(r_ECEF_actual);
            disp(r_ECEF_expected);
            disp('ECEF Position Comparison (Expected vs Reference):');
            disp(r_ECEF_expected);
            disp(r_ECEF_ref);
            disp('ECEF Velocity Comparison (Actual vs Reference):');
            disp(v_ECEF_actual);
            disp(v_ECEF_ref);
            
            testCase.verifyEqual(r_ECI_actual, r_ECI_builtin, 'AbsTol',.30); % 30 cm
            testCase.verifyEqual(r_ECI_builtin, r_ECI_ref, 'AbsTol',.30);
            testCase.verifyEqual(v_ECI_actual,v_ECI_ref,'AbsTol', .30); % 30 cm/sec
            testCase.verifyEqual(r_ECEF_actual, r_ECEF_expected, 'AbsTol',.30); % 30 cm
            testCase.verifyEqual(r_ECEF_expected, r_ECEF_ref, 'AbsTol',.30);
            testCase.verifyEqual(v_ECEF_actual, v_ECEF_ref, 'AbsTol', .30); % 30 cm/sec
        end

        function test_ECI_to_ECEF_Project(testCase)
            [S,ENV] = Setup.loadSettings('F','Final_3D',false,false);


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
            
            disp('ECEF Vectors'); disp(r_ECEF_actual); disp(r_ECEF_expected);
            disp('rot matrices diff')
            disp(R_ECEF_from_ECI_actual - R_ECI_to_ECEF);
            
            % diff_v_m_per_s = abs(v_ECEF_actual - v_ECEF_expected)

            testCase.verifyEqual(r_ECEF_actual,r_ECEF_expected,'AbsTol',1.5);
            testCase.verifyEqual(R_ECEF_from_ECI_actual, R_ECI_to_ECEF,'AbsTol',1e-4);
            % testCase.verifyLessThanOrEqual(diff_v_m_per_s, 1e-1);
        end
    end
end