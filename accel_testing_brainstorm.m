clear; clc;
%%
case_names = {'A','B','C','D','E','F'};
scenario_config = {'Accel','Final_3D','Final_1D','HW5'};
[S,ENV] = Setup.loadSettings('F','Accel',true,true);

%%

STM_0 = eye(7);

X_state = [S.IC_Sat_Epoch.position_ECI_meters;...
           S.IC_Sat_Epoch.velocity_ECI_meters_per_second;...
           1.88;...
           STM_0(:)...
           ];
time_epoch_struct = Tools.ComputeTimeSystems(S.IC_Sat_Epoch.epoch_date_time_UTC);

%%

[r0_ECI_station1_km, v0_ECI_station1_km_s] = ECEF_ECI_Converter(S.S, v0_ECEF_station_1_km, S.IC_Sat_Epoch.epoch_date_time_UTC, "ECEF_to_ECI", EOP_interpolated);
[r0_ECI_station2_km, v0_ECI_station2_km_s] = ECEF_ECI_Converter(r0_ECEF_station_2_km, v0_ECEF_station_2_km, S.IC_Sat_Epoch.epoch_date_time_UTC, "ECEF_to_ECI", EOP_interpolated);
[r0_ECI_station3_km, v0_ECI_station3_km_s] = ECEF_ECI_Converter(r0_ECEF_station_3_km, v0_ECEF_station_3_km, S.IC_Sat_Epoch.epoch_date_time_UTC, "ECEF_to_ECI", EOP_interpolated);

%%
clc;
dX_dt = jah_sat_1_ode(0,X_state,S, ENV,true);

%% Time how long this ode takes to call
tic;
for i = 1:1000
    jah_sat_1_ode(0,X_state, time_epoch_struct.jd_UTC_days,false);
end

a = toc/1000;
disp(['Average time per call: ', num2str(a), ' seconds']);

