clear; clc;
%%
[S,ENV] = Setup.loadSettings('F','1 Hour','Final - 1 Day Subset',false,true);

STM_0 = eye(7);
X_state_epoch = [S.IC_Sat_Epoch.position_ECI_meters;...
                   S.IC_Sat_Epoch.velocity_ECI_meters_per_second;...
                   1.88;...
                   STM_0(:)...
                   ];
time_epoch_struct = Tools.ComputeTimeSystems(S.IC_Sat_Epoch.epoch_date_time_UTC);
ENV.EOP_IERS.MJD_days = time_epoch_struct.mjd_UTC_days;
options = odeset('RelTol',3e-10,'AbsTol',1e-12); % TODO: Look at these more closely.
%%
% A priori information
P_bar_cov_0 = S.StateCovariances;

t_obs = S.ref_data.Actual_Measurements.time_sec_past_epoch;
N_obs = length(t_obs);

%% Create storage vector for computed observations!

Y_prefit_computed_observations = zeros(2,N_obs); % To store computed observations for each time step
Y_postfit_computed_observations = zeros(2,N_obs);

%% Initialize i_minus_1 (ONLY FOR FIRST TIME!)
t_i_minus_1 = t_obs(1);
X_star_input_ti_minus_1 = X_state_epoch;
P_cov_i_minus_1 = P_bar_cov_0;

%% Initial Observation at t=0
if t_obs(1) == 0 
    % 1. Standard Setup (Always do this)
    station_id = ENV.ref_data.Actual_Measurements.station_id(1);
    Y_1 = ENV.ref_data.Actual_Measurements{1, {'apparent_range_km', 'apparent_range_rate_km_s'}}(:);
    R_1 = IC.Stations(station_id).Covariance;
    
    date_time_0 = S.IC_Sat_Epoch.epoch_date_time_UTC;
    EOP_params_0 = Tools.interpolate_EOP(date_time_0, ENV.EOP_IERS, ENV.EOP_Celestrak);
    [r_stn_t_0_ECI_meters, v_stn_t_0_ECI_meters, ~, ~] = Tools.ECEF_ECI_Converter(S.Stations(station_id).position_ECEF_meters, zeros(3,1), date_time_0, "ECEF_to_ECI", EOP_params_0);
    % 2. Pre-fit Residual (Always based on the original IC)
    Y_comp_0 = Measurements.Compute_Range_Range_Rate(X_star_0(1:3), X_star_0(4:6), r_stn_t_0_ECI_meters, v_stn_t_0_ECI_meters);
    Y_prefit_computed_observations(:,1) = Y_comp_0;
    y_residuals_0 = Y_1 - Y_comp_0;

    % 3. Check if this station is "Active"
    update_t0 = false;
    if strcmp(active_station_mode, 'All') || (station_id == station_map.(active_station_mode))
        update_t0 = true;
    end

    if update_t0
        % --- UPDATE MODE ---
        H_0 = ENV.HmatrixFcn(X_star_0(1:3), X_star_0(4:6), r_stn_t_0_ECI_meters, v_stn_t_0_ECI_meters);
        K_0 = P_bar_cov_0 * H_0' / (H_0 * P_bar_cov_0 * H_0' + R_1);
        
        % Update State and Covariance
        x_hat_0 = K_0 * y_residuals_0; 
        X_corrected_0 = X_star_0 + x_hat_0;
        I = eye(N_states);
        P_cov_i_minus_1 = (I - K_0 * H_0) * P_bar_cov_0 * (I - K_0 * H_0)' + K_0 * R_1 * K_0';
    else
        % --- PASSIVE MODE ---
        % No correction applied to the state
        X_corrected_0 = X_star_0; 
        % Covariance remains the A-Priori P
        P_cov_i_minus_1 = P_bar_cov_0; 
    end

    % 4. Post-fit calculation (Using whatever the state is now)
    Y_post = Measurements.Compute_Range_Range_Rate(X_corrected_0(1:3), X_corrected_0(4:6), r_stn_t_0_ECI_meters, v_stn_t_0_ECI_meters);
    Y_postfit_computed_observations(:, 1) = Y_post;
    
    % 5. Prepare for the first ODE step
    X_star_input_ti_minus_1 = [X_corrected_0; STM_t0_t0(:)];
    
    fprintf('t=0 (i=1) Processing Complete. Update Applied: %s\n', string(update_t0));
end

%% Start for-loop
for i = 2:length(t_obs)

    % CHECK IF WE SHOULD UPDATE THIS STEP!
    update_this_step = false;
    if strcmp(active_station_mode, 'All')
        update_this_step = true;
    elseif station_id == station_map.(active_station_mode)
        update_this_step = true;
    end

%% Read Next Observation Y_i and R_i

station_id = S.ref_data.Actual_Measurements.station_id(i);
actual_obs_y_i = S.ref_data.Actual_Measurements{i, {'apparent_range_km', 'apparent_range_rate_km_s'}}(:);
R_conv_i = S.Stations(i).Covariance;

%% Propogate to Observation Time t_1

[t,y] = ode45(@(t,X) jah_sat_1_ode(t, X, time_epoch_struct.jd_UTC_days, false), ...
            [t_obs(i-1),t_obs(i)], X_star_input_ti_minus_1, options);

r_sat_t_i_ECI_meters = y(end,1:3)';
v_sat_t_i_ECI_meters_s = y(end,4:6)';
full_orbit_ECI = y(:,1:6);          % For plotting the orbit path in ECI
X_star_ti = y(end,1:N_states)';
STM_t_i_t_i_minus_1 = reshape(y(end,N_states+1:end), N_states, N_states);

% Visuals.plot_position(r_sat_t_i_ECI_meters, v_sat_t_i_ECI_meters_s, P_cov_i_minus_1, full_orbit_ECI);

%% Time Update of State Estimate and Covariance + PROCESS NOISE!
% x_bar_correction_i = STM_t_i_t_i_minus_1 * x_hat_correction_i_minus_1; %
% Omitted from EKF

State_Noise_Compensation_7x7 = zeros(N_states, N_states); % Pad to 7x7 (assuming no process noise on the 7th state, Cd)

% delta_t = t_obs(i) - t_obs(i-1);
% Q_matrix = eye(3) * IC.Sigma_Accel_meters_s2^2;
% State_Noise_Compensation_6x6 = delta_t^2 *  [(delta_t^2/4)*Q_matrix, (delta_t/2)*Q_matrix;
%                                              (delta_t/2)*Q_matrix, (1)*Q_matrix ];

% State_Noise_Compensation_7x7(1:6, 1:6) = State_Noise_Compensation_6x6; %
% TODO -- Turn process noise model back on later.

P_bar_cov_i = STM_t_i_t_i_minus_1 * P_cov_i_minus_1 * STM_t_i_t_i_minus_1' + State_Noise_Compensation_7x7;

%% Compute observation deviation, observation state matrix, gain matrix

date_time_i =  time_epoch_struct.jd_UTC_days + seconds(t_obs(i));
EOP_params_i = Tools.interpolate_EOP(date_time_i, ENV.EOP_IERS, ENV.EOP_Celestrak);
[r_station_t_i_ECI_meters, v_station_t_i_ECI_meters] = Tools.ECEF_ECI_Converter(S.Stations(station_id).position_ECEF_meters,zeros(3,1),date_time_i,"ECEF_to_ECI",EOP_params_i);
computed_obs_y_i = Measurement.Compute_Range_Range_Rate(r_sat_t_i_ECI_meters,v_sat_t_i_ECI_meters_s,r_station_t_i_ECI_meters,v_station_t_i_ECI_meters);
Y_prefit_computed_observations(:,i) = computed_obs_y_i;
y_residuals_i = actual_obs_y_i - computed_obs_y_i;

%% Kalman Update
if update_this_step
    H_tilde_i = computed_obs_y_i;
    K_i = P_bar_cov_i * H_tilde_i' / (H_tilde_i * P_bar_cov_i * H_tilde_i' + R_conv_i);
    x_hat_correction_i = K_i * (y_residuals_i); % This is EKF unique
    X_star_ti = X_star_ti + x_hat_correction_i;

    I = eye(N_states);
    % P_cov_i = (I - K_i * H_tilde_i) * P_bar_cov_i;
    % Use the Joseph Form to ensure PSD!
    % Measurement Update of Covariance (Joseph Form)
    P_cov_i = (I - K_i * H_tilde_i) * P_bar_cov_i * (I - K_i * H_tilde_i)' + K_i * R_conv_i * K_i';
    P_cov_i = 0.5 * (P_cov_i + P_cov_i');
else
    P_cov_i = P_bar_cov_i; % No update to covariance if not updating this step
    x_hat_correction_i = zeros(N_states,1); % No correction if not updating
end


%% Post-fit calculation
% 1. Get the updated total state
% In your EKF, X_star_ti has already been updated with the correction!
X_updated = X_star_ti; 

% 2. Re-compute observation with updated state
Y_post = Measurements.Compute_Range_Range_Rate(X_updated(1:3), X_updated(4:6), r_station_t_i_ECI_meters, v_station_t_i_ECI_meters);
Y_postfit_computed_observations(:, i) = Y_post;


%% Print Results
if print_updates || (mod(i,20) == 0)
    disp('State Estimate Correction at t_i:');
    disp(x_hat_correction_i);
    disp('Covariance Matrix at t_i:');
    disp(P_cov_i);
    fprintf('Step %d of %d complete. Time: %.2f sec\n', i, length(t_obs), t_i);
    disp('State Correction:'); disp(x_hat_correction_i');
    
    % This will pause execution until you press any key in the Command Window
    % disp('Press any key to continue to the next observation...');
    % pause;

end

%% Update for next iteration

% Extract the nominal state from the end of the current ODE integration
STM_reset = eye(N_states);
X_star_input_ti_minus_1 = [X_star_ti;STM_reset(:)];

% Pass the updated covariance and state deviation to the next step
P_cov_i_minus_1 = P_cov_i;
% x_hat_correction_i_minus_1 = x_hat_correction_i; % Shouldn't need this
% actually

%% END FOR LOOP
end

