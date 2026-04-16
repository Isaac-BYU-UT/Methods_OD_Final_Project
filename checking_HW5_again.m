clear; clc;

%% SETUP
case_names = {'A','B','C','D','E','F'};
scenario_config = {'Accel','Final_3D','Final_1D','HW5'};
[S, ENV] = Setup.loadSettings('F','HW5',false,false);

%%

initial_state = EKF.initialize_ekf(S, ENV);

options = odeset('RelTol',1e-9,'AbsTol',1e-11);
[t, y] = ode45(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, false), ...
        [0, 21600], initial_state.X_input, options);
%%
(reshape(y(end,8:56),7,7) - S.ref_data.STM_21600_0_ref)./S.ref_data.STM_21600_0_ref