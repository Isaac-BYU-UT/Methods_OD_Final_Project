clear; clc;
%%
[S,ENV] = Setup.loadSettings('F','1 Hour','Acceleration Testing A');

%%

STM_0 = eye(7);

X_state = [S.IC_Sat_Epoch.position_ECI_meters;...
           S.IC_Sat_Epoch.velocity_ECI_meters_per_second;...
           1.88;...
           STM_0(:)...
           ];
time_epoch_struct = Tools.ComputeTimeSystems(S.IC_Sat_Epoch.epoch_date_time_UTC);

%%
clc;
dX_dt = jah_sat_1_ode(0,X_state,time_epoch_struct.jd_UTC_days,true);

%% Time how long this ode takes to call
tic;
for i = 1:1000
    jah_sat_1_ode(0,X_state, time_epoch_struct.jd_UTC_days,false);
end

a = toc/1000;
disp(['Average time per call: ', num2str(a), ' seconds']);

