clear; clc;

%% SETUP
i_case = 6;
case_names = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

i_scenario = 2;
scenario_config = {'Accel', 'Final_1D', 'Final_3D', 'Final_6D', 'HW5', '24Dynamics'};

[S, ENV] = Setup.loadSettings(case_names{i_case}, scenario_config{i_scenario}, false, false);
options = odeset('RelTol',1e-12,'AbsTol',1e-14);
%% Initialize Filter Values

X_star_0 = [S.IC_Sat_Epoch.position_ECI_meters;S.IC_Sat_Epoch.velocity_ECI_meters_per_second];
N_states = length(X_star_0);
STM_t0_t0 = eye(N_states);
y0 = [X_star_0; STM_t0_t0(:)]; % For input into ODE

states_initial_input = y0;
t_obs = S.ref_data.Actual_Measurements.time_sec_past_epoch;

t_obs = t_obs(t_obs < 6 * Units.SEC_IN_SOLAR_DAY/24);
threshold_position_m = 1; % Stops changing by max 100 meters;
state_change_position_m = 1e10; % initialize high

threshold_velocity_m_s = .1; % m/s
state_change_velocity_m_s = 1e10; % initialize high

iteration = 0;
while (state_change_position_m > threshold_position_m || state_change_velocity_m_s > threshold_velocity_m_s)
%% Propagation of Reference Trajectory
disp(iteration);
iteration = iteration + 1;
disp(states_initial_input(1:N_states));

[t,y] = ode45(@(t,X) jah_sat_1_ode( ...
    t, X, S, ENV, false), ...
    t_obs, states_initial_input, options);


%% Extract States
% Extract r_sat and v_sat at the observation times from the ODE solution
X_star_table = array2table(y(:,1:N_states), 'VariableNames', {'r_x_m', 'r_y_m', 'r_z_m', 'v_x_m_s', 'v_y_m_s', 'v_z_m_s'}); % Convert to table for easier handling
X_star_table.time_sec_past_epoch = t; % Add time column to the table

% Now add STM elements to table.
X_star_table.STM_elements = cell(height(X_star_table),1); % Pre-allocate cell array for STM elements
for i = 1:height(X_star_table)
    STM_i = reshape(y(i,N_states+1:end), N_states, N_states); % Extract STM elements for this time step
    X_star_table.STM_elements{i} = STM_i; % Store in table
end

%% H Matrix and STM at Each Time Step
N_obs = length(t_obs);
STM_i = NaN(N_states, N_states, N_obs); % Pre-allocate 3D array for STM at each time step
H_i = NaN(2, N_states, N_obs); % Pre-allocate 3D array for H at each time step
H_mapped_epoch = NaN(2, N_states, N_obs);
residual_prefit = NaN(2, N_obs); % 2 x N_obs to make column vectors more naturally
Y_prefit_computed = NaN(2, N_obs); % 2 x N_obs to make column vectors more naturally

% A priori information
P0 = S.StateCovariances.P_Covariance_States(1:N_states, 1:N_states); % a priori covariance
W_apriori = P0 \ eye(N_states); % This is W_bar_k

% Initialized Accumulated Normal Equations
% Information_Matrix = zeros(N_states, N_states);
Information_Matrix = W_apriori;
Residual_Vector = zeros(N_states, 1); % The apriori estimate for x_0 is just the 0 vector

for i = 1:length(t_obs)
    
    STM_i(:,:,i) = X_star_table.STM_elements{i}; % Extract STM for this time step and store in 3D array

    % Station, Time, EOP_parameters for table row
    station_id = S.ref_data.Actual_Measurements.station_id(i);
    date_time_i = S.IC_Sat_Epoch.epoch_date_time_UTC + seconds(t_obs(i));
    EOP_params_i = Tools.interpolate_EOP(date_time_i, ENV.EOP_IERS, ENV.EOP_Celestrak);

    [r_station_t_i_ECI_m, v_station_t_i_ECI_m] = Tools.ECEF_ECI_Converter(S.Stations(station_id).position_ECEF_meters, ...
                                                                                zeros(3,1), date_time_i, "ECEF_to_ECI", EOP_params_i);

    
    H_i(:,:,i) = Measurements.Compute_H_matrix( ...
        X_star_table{i,1:3}(:), ... % r_sat_ECI_m
        X_star_table{i,4:6}(:), ... % v_sat_ECI_m_s 
        r_station_t_i_ECI_m(:), ...
        v_station_t_i_ECI_m(:));

    H_mapped_epoch(:,:,i) = H_i(:,:,i)*STM_i(:,:,i); % dG/dX |t_i * STM(t_i, t_0) to map from epoch to time of measurement
    
    % Compute measurement at time step i
    Y_prefit_computed(:,i) = Measurements.Compute_Range_Range_Rate( ...
        X_star_table{i,1:3}(:), ... % r_sat_ECI_m
        X_star_table{i,4:6}(:), ... % v_sat_ECI_m_s 
        r_station_t_i_ECI_m(:), ...
        v_station_t_i_ECI_m(:))';

    y_obs_meters_i =  [S.ref_data.Actual_Measurements.apparent_range_km(i); ...
                               S.ref_data.Actual_Measurements.apparent_range_rate_km_s(i)] ...
                                * Units.KILOMETERS; % SUPER IMPORTANT!!!

    % Compute REsidual
    residual_prefit(:,i) = y_obs_meters_i - Y_prefit_computed(:,i); % THESE ARE COLUMN VECTORS!

    % Note that the \eye(2) is to compute the inverse of the covariance matrix:
    R_cov_inv = S.Stations(station_id).Covariance\eye(2); % Inverse of covariance matrix for this station
    Information_Matrix = Information_Matrix + H_mapped_epoch(:,:,i)' * R_cov_inv * H_mapped_epoch(:,:,i); % Accumulate Information Matrix
    Information_Matrix = 0.5 * (Information_Matrix + Information_Matrix'); % INSURE THAT INFORMATION MATRIX REMAINS SYMMETRIC!
    Residual_Vector = Residual_Vector + H_mapped_epoch(:,:,i)' * R_cov_inv * residual_prefit(:,i); % Accumulate Residual Vector
end

%% Solve at Epoch
x_correction_0 = Information_Matrix\Residual_Vector; % Solve for correction at epoch
X_star_new_0 = X_star_0 + x_correction_0; % Update initial state estimate at epoch

P_estimated = Information_Matrix \ eye(N_states); % Updated estimate of state covariance at epoch (inverse of information matrix)

%% Post Fit Residuals
epsilon_postfit_i = NaN(2, N_obs); % 2 x N_obs to make column vectors more naturally
Y_postfit_computed = NaN(2,N_obs);
for i = 1:length(t_obs)
    epsilon_postfit_i(:,i) = residual_prefit(:,i) - H_mapped_epoch(:,:,i)*x_correction_0; % Post-fit residuals at time step i
    Y_postfit_computed(:,i) = Y_prefit_computed(:,i) + H_mapped_epoch(:,:,i)*x_correction_0;
end

%% Plot Pre-fit Residuals
% 
% Measurement_Table = Visuals.make_measurement_table(t_obs, S.ref_data.Actual_Measurements, Y_prefit_computed');
% Visuals.plot_station_residuals_no_bounds(Measurement_Table, {S.Stations.name});
% 
% %% Now post-fit residuals
% 
% Measurement_Table = Visuals.make_measurement_table(t_obs, S.ref_data.Actual_Measurements, Y_postfit_computed');
% Visuals.plot_station_residuals_no_bounds(Measurement_Table, {S.Stations.name});

%% Update ODE Run
X_star_0 = X_star_new_0;
states_initial_input = [X_star_new_0; STM_t0_t0(:)];
state_change_position_m = max(abs(x_correction_0(1:3)));
fprintf("Max Position Change: %d m", state_change_position_m);

state_change_velocity_m_s = max(abs(x_correction_0(4:6)));
fprintf("Max Vel Change: %d m/s", state_change_velocity_m_s);
end

