% CONVENTIONAL SEQUENTIAL KALMAN FILTER
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
N_obs = length(t_obs);

%% Create storage vector for computed observations!

Y_prefit_computed_observations = zeros(2,N_obs); % To store computed observations for each time step
Y_postfit_computed_observations = zeros(2,N_obs);

%% Initialize i_minus_1 (ONLY FOR FIRST TIME!)
t_i_minus_1 = t_obs(1);
X_star_input_ti_minus_1 = IC.y0;
P_cov_i_minus_1 = P_bar_cov_0;
x_hat_correction_i_minus_1 = x_bar_correction_0;


%% Initial Observation at t=0 (First Step)
% Check if the first observation is indeed at the epoch
if t_obs(1) == 0  % Or whatever your epoch time is
    %% Conventional
    % Initial Observation at t=0 (Conventional/Linearized KF)
    % In a conventional/linearized filter, we do NOT update X_star. 
    % We only update the state deviation (x_hat) and covariance (P).

    % 1. Get observation data and station info for t=0
    station_id = ENV.ref_data.Actual_Measurements.station_id(1);
    Y_1 = ENV.ref_data.Actual_Measurements{1, {'apparent_range_km', 'apparent_range_rate_km_s'}}(:);
    R_1 = IC.Stations(station_id).Covariance;
    
    date_time_0 = IC.UTC_epoch;
    EOP_params_0 = Tools.interpolate_EOP(date_time_0, "IERS");
    [r_stn_0, v_stn_0] = Tools.ECEF_ECI_Converter(IC.Stations(station_id).r_ECEF_km, ...
                         zeros(3,1), date_time_0, "ECEF_to_ECI", EOP_params_0);

    % 2. Compute the observation matrix (H) and prefit residual (y)
    % We use the initial reference state (X_star_0)
    Y_comp_0 = ENV.MeasFcn(X_star_0(1:3), X_star_0(4:6), r_stn_0, v_stn_0);
    H_0 = ENV.HmatrixFcn(X_star_0(1:3), X_star_0(4:6), r_stn_0, v_stn_0);
    
    % Prefit residual (Innovation)
    y_residuals_0 = Y_1 - Y_comp_0;
    Y_prefit_computed_observations(:,1) = Y_comp_0;

    % 3. Kalman Gain
    K_0 = P_bar_cov_0 * H_0' / (H_0 * P_bar_cov_0 * H_0' + R_1);

    % 4. Update State Correction and Covariance
    % Formula: x_hat = x_bar + K * (y - H*x_bar)
    % At t=0, x_bar (x_bar_correction_0) is typically zeros.
    x_hat_correction_i_minus_1 = x_bar_correction_0 + K_0 * (y_residuals_0 - H_0 * x_bar_correction_0);

    % 5. Update Covariance (Joseph Form)
    I = eye(N_states);
    P_cov_i_minus_1 = (I - K_0 * H_0) * P_bar_cov_0 * (I - K_0 * H_0)' + K_0 * R_1 * K_0';
    
    % CRITICAL FOR CONVENTIONAL/LKF: 
    % We do NOT add x_hat to X_star. The ODE will continue to use the original 
    % uncorrected X_star_0 as the integration base.
    X_star_input_ti_minus_1 = [X_star_0; STM_t0_t0(:)]; 

    %% Post-fit calculation
    % 1. Get the updated total state
    % Note: for EKF this is X_star_ti (already updated), for CKF it is X_star_ti + x_hat
    X_updated = X_star_0 + x_hat_correction_i_minus_1; 
    
    % 2. Re-compute observation with updated state
    % We keep the station position/velocity from the pre-fit calculation
    Y_post = ENV.MeasFcn(X_updated(1:3), X_updated(4:6), r_stn_0, v_stn_0);
    Y_postfit_computed_observations(:, 1) = Y_post; % Use index 1, not i
    
    fprintf('t=0 (i=1) Conventional Update Complete.\n');
end

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
computed_observation_i = ENV.MeasFcn(r_sat_t_i_ECI_km,v_sat_t_i_ECI_km_s,r_station_t_i_ECI_km,v_station_t_i_ECI_km);
Y_prefit_computed_observations(:,i) = computed_observation_i;
y_residuals_i = Y_i - computed_observation_i;
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

% Extract sigma for Range and Range-rate from the diagonal of R
sigma_range = sqrt(R_conv_i(1,1));
sigma_range_rate = sqrt(R_conv_i(2,2));

% Store 3-sigma values for plotting
three_sigma_bounds(:, i) = [3 * sigma_range; 3 * sigma_range_rate];

%% Post-fit calculation
% 1. Get the updated total state
% Note: for EKF this is X_star_ti (already updated), for CKF it is X_star_ti + x_hat
X_updated = X_star_ti + x_hat_correction_i; 

% 2. Re-compute observation with updated state
% We keep the station position/velocity from the pre-fit calculation
Y_post = ENV.MeasFcn(X_updated(1:3), X_updated(4:6), r_station_t_i_ECI_km, v_station_t_i_ECI_km);
Y_postfit_computed_observations(:, i) = Y_post;

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

%% Plot Position RIC graphs
Visuals.plot_position(r_sat_t_i_ECI_km, v_sat_t_i_ECI_km_s, P_cov_i_minus_1, full_orbit_ECI);

%% Plot Pre-fit Residuals (Y - h(X_bar)) --> How did dynamic model do?

Measurement_Table = Visuals.make_measurement_table(t_obs,ENV.ref_data.Actual_Measurements,Y_prefit_computed_observations'); % Ensure correct dimensions!
Visuals.plot_station_residuals(Measurement_Table, {IC.Stations.Name});
Visuals.plot_measurement_correlation_linked(Measurement_Table);

%% Plot Post-fit Residuals (Y - h(X_hat)) --> How well did filter do?
Post_Measurement_Table = Visuals.make_measurement_table(t_obs, ENV.ref_data.Actual_Measurements, Y_postfit_computed_observations');
Visuals.plot_station_residuals_3_sigma(Post_Measurement_Table, {IC.Stations.Name}, IC);
Visuals.plot_measurement_correlation_linked(Measurement_Table);

%% Compute RMS of pre-fit and post-fit residuals
range_residuals = Measurement_Table.apparent_range_km - Measurement_Table.computed_range_km;
range_rate_residuals = Measurement_Table.apparent_range_rate_km_s - Measurement_Table.computed_range_rate_km_s;
RMS_range_prefit = sqrt(mean(range_residuals.^2));
RMS_range_rate_prefit = sqrt(mean(range_rate_residuals.^2));
range_residuals_post = Post_Measurement_Table.apparent_range_km - Post_Measurement_Table.computed_range_km;
range_rate_residuals_post = Post_Measurement_Table.apparent_range_rate_km_s - Post_Measurement_Table.computed_range_rate_km_s;
RMS_range_postfit = sqrt(mean(range_residuals_post.^2));  
RMS_range_rate_postfit = sqrt(mean(range_rate_residuals_post.^2));  

fprintf('RMS of Range Residuals (Prefit): %.6f km\n', RMS_range_prefit);
fprintf('RMS of Range Rate Residuals (Prefit): %.6f km/s\n', RMS_range_rate_prefit);
fprintf('RMS of Range Residuals (Postfit): %.6f km\n', RMS_range_postfit);
fprintf('RMS of Range Rate Residuals (Postfit): %.6f km/s\n', RMS_range_rate_postfit);


