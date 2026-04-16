[S,ENV] = Setup.loadSettings('F','1 Hour','Final - Full 3 Days');

tic;
for i = 1:100

[r_ECEF_actual,v_ECEF_actual,~,~] = Tools.ECEF_ECI_Converter(S.IC_Sat_Epoch.position_ECI_meters,...
                                                            S.IC_Sat_Epoch.velocity_ECI_meters_per_second, ...
                                                            S.IC_Sat_Epoch.epoch_date_time_UTC,...
                                                            'ECEF_to_ECI',...
                                                            ENV.EOP_t0);

end

t1 = toc;

for i = 1:100
    R_ECI_to_ECEF = dcmeci2ecef('IAU-76/FK5',...
                                S.IC_Sat_Epoch.epoch_date_time_UTC,...
                                ENV.EOP_t0.delta_AT_sec, ...
                                ENV.EOP_t0.UT1_UTC_sec,...
                                [ENV.EOP_t0.x_pole_arcsec, ENV.EOP_t0.y_pole_arcsec] * Units.ARCSEC_TO_RAD,...
                                'dNutation',[ENV.EOP_t0.dPsi_milli_arcsec,ENV.EOP_t0.dEpsilon_milli_arcsec] * Units.MILLI_TO_NOM * Units.ARCSEC_TO_RAD...
                                );
end

t2 = toc - t1;

for i = 1:100

    r_ECEF_expected = eci2ecef(S.IC_Sat_Epoch.epoch_date_time_UTC,...
                                S.IC_Sat_Epoch.position_ECI_meters,...
                                S.IC_Sat_Epoch.velocity_ECI_meters_per_second,...
                                [0;0;0], ...
                                'dAT',ENV.EOP_t0.delta_AT_sec,...
                                'dUT1', ENV.EOP_t0.UT1_UTC_sec,...
                                'pm',[ENV.EOP_t0.x_pole_arcsec, ENV.EOP_t0.y_pole_arcsec] * Units.ARCSEC_TO_DEG ...
                                );

end

t3 = toc - t2;

disp([t1, t2, t3])

% 
% xp_rad = ENV.EOP_t0.x_pole_arcsec * Units.ARCSEC_TO_RAD;
% yp_rad = ENV.EOP_t0.y_pole_arcsec * Units.ARCSEC_TO_RAD;
% 
% % Approximate Version
% M = [1 0  -xp_rad;
%      0 1    yp_rad;
%      xp_rad -yp_rad 1];

r_ECEF_expected = R_ECI_to_ECEF \ S.IC_Sat_Epoch.position_ECI_meters;

% v_ECEF_expected = 0;

% Omega_Earth_Vector_Rad_Sec = [0; 0; EOP_params.omega_earth_rad_sec];
% TODO: Come back to this later!

diff_r_meters = abs(r_ECEF_actual - r_ECEF_expected)
% diff_v_m_per_s = abs(v_ECEF_actual - r_ECEF_expected)

% testCase.verifyLessThanOrEqual(diff_r_meters,1e-1);
% testCase.verifyLessThanOrEqual(diff_v_m_per_s, 1e-1);