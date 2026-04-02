clear; clc;
% This is just section 4.7 from Schutz "Statistical Orbit Determination"!
%% Load Initial Conditions
IC = Initial_Conditions.initial_conditions_HW5();
%% Setup Environment
ENV = Environment.setup_environment(IC);
options = odeset('RelTol',3e-14,'AbsTol',1e-16);
%% Initialize Filter Values
X_star_0 = [IC.r_sat_ECI_km; IC.v_sat_ECI_km_s; IC.C_drag];
N_states = length(X_star_0);
STM_t0_t0 = eye(N_states);
IC.y0 = [X_star_0; STM_t0_t0(:)]; % For input into ODE

% A priori information
P_bar_cov_0 = IC.P_Covariance_States; % This is the covariance matrix of the initial state estimate at epoch (7x7)
% TODO: ADD STATE PROCESS NOISE!!!
x_bar_correction_0 = zeros(N_states,1);
%% Set Initial Observation Times (ONLY FOR FIRST ROUND!)
t_obs = ENV.ref_data.Actual_Measurements.time_sec_past_epoch;

%% Initialize i_minus_1 (ONLY FOR FIRST TIME!)
t_i_minus_1 = t_obs(1);
X_star_input_ti_minus_1 = IC.y0;
P_cov_i_minus_1 = P_bar_cov_0;
x_hat_correction_i_minus_1 = x_bar_correction_0;

%% Start for-loop
for i = 2:length(t_obs)

%% Read Next Observation Y_i and R_i

t_i = t_obs(i);
t_i_minus_1 = t_obs(i-1);
station_id = ENV.ref_data.Actual_Measurements.station_id(i);
Y_i = ENV.ref_data.Actual_Measurements{i, {'apparent_range_km', 'apparent_range_rate_km_s'}}(:);
R_conv_i = IC.Stations(station_id).Covariance;

%% Propogate to Observation Time t_1

[t,y] = ode45(@(t,X) jah_sat_1_ode( ...
    t, X, ENV.AccelFcn, ENV.AmatrixFcn, IC.time_struct.jd_UTC_days), ...
    [t_i_minus_1,t_i], X_star_input_ti_minus_1, options);

r_sat_t_i_ECI_km = y(end,1:3)';
v_sat_t_i_ECI_km_s = y(end,4:6)';
full_orbit_ECI = y(:,1:6); % For plotting the orbit path in ECI
X_star_ti = y(end,1:N_states)';
STM_t_i_t_i_minus_1 = reshape(y(end,N_states+1:end), N_states, N_states);

% Visuals.plot_position(r_sat_t_i_ECI_km, v_sat_t_i_ECI_km_s, P_cov_i_minus_1, full_orbit_ECI);

%% Time Update of State Estimate and Covariance + PROCESS NOISE!
x_bar_correction_i = STM_t_i_t_i_minus_1 * x_hat_correction_i_minus_1;

delta_t = t_i - t_i_minus_1;
Q_matrix = eye(3) * IC.Sigma_Accel_km_s2^2;

State_Noise_Compensation_6x6 = delta_t^2 *  [(delta_t^2/4)*Q_matrix, (delta_t/2)*Q_matrix;
                                             (delta_t/2)*Q_matrix, (1)*Q_matrix ];

State_Noise_Compensation_7x7 = zeros(N_states, N_states); % Pad to 7x7 (assuming no process noise on the 7th state, Cd)
State_Noise_Compensation_7x7(1:6, 1:6) = State_Noise_Compensation_6x6;


P_bar_cov_i = STM_t_i_t_i_minus_1 * P_cov_i_minus_1 * STM_t_i_t_i_minus_1' + State_Noise_Compensation_7x7;

%% Compute observation deviation, observation state matrix, gain matrix

date_time_i = IC.UTC_epoch + seconds(t_i);
EOP_params_i = Tools.interpolate_EOP(date_time_i, "IERS");
[r_station_t_i_ECI_km, v_station_t_i_ECI_km] = Tools.ECEF_ECI_Converter(IC.Stations(station_id).r_ECEF_km,zeros(3,1),date_time_i,"ECEF_to_ECI",EOP_params_i);
y_residuals_i = Y_i - ENV.MeasFcn(r_sat_t_i_ECI_km,v_sat_t_i_ECI_km_s,r_station_t_i_ECI_km,v_station_t_i_ECI_km);
H_tilde_i = ENV.HmatrixFcn(r_sat_t_i_ECI_km,v_sat_t_i_ECI_km_s,r_station_t_i_ECI_km,v_station_t_i_ECI_km);
K_i = P_bar_cov_i * H_tilde_i' / (H_tilde_i * P_bar_cov_i * H_tilde_i' + R_conv_i);

%% Measurement Update of State Estimate and Covariance
x_hat_correction_i = x_bar_correction_i + K_i * (y_residuals_i - H_tilde_i * x_bar_correction_i);
I = eye(N_states);
% P_cov_i = (I - K_i * H_tilde_i) * P_bar_cov_i;
% Use the Joseph Form to ensure PSD!
% Measurement Update of Covariance (Joseph Form)
P_cov_i = (I - K_i * H_tilde_i) * P_bar_cov_i * (I - K_i * H_tilde_i)' + K_i * R_conv_i * K_i';
P_cov_i = 0.5 * (P_cov_i + P_cov_i');


%% Print Results
disp('State Estimate Correction at t_i:');
disp(x_hat_correction_i);
disp('Covariance Matrix at t_i:');
disp(P_cov_i);

%% Print and Plot Results
fprintf('Step %d of %d complete. Time: %.2f sec\n', i, length(t_obs), t_i);
disp('State Correction:'); disp(x_hat_correction_i');

% This will pause execution until you press any key in the Command Window
% disp('Press any key to continue to the next observation...');
% pause;

%% Update for next iteration

% Extract the nominal state from the end of the current ODE integration
STM_reset = eye(N_states);
X_star_input_ti_minus_1 = [X_star_ti;STM_reset(:)];

% Pass the updated covariance and state deviation to the next step
P_cov_i_minus_1 = P_cov_i;
x_hat_correction_i_minus_1 = x_hat_correction_i;

%% END FOR LOOP
end

%% Plots and Analysis
Visuals.plot_position(r_sat_t_i_ECI_km, v_sat_t_i_ECI_km_s, P_cov_i_minus_1, full_orbit_ECI);

