clear; clc;
%% Load Initial Conditions
IC = Initial_Conditions.initial_conditions_HW5();

%% Setup Environment
ENV = Environment.setup_environment(IC);
%% Initialize Filter Values

X_star_0 = [IC.r_sat_ECI_km; IC.v_sat_ECI_km_s; IC.C_drag];
N_states = length(X_star_0);
STM_t0_t0 = eye(N_states);
IC.y0 = [X_star_0; STM_t0_t0(:)]; % For input into ODE

states_initial_input = IC.y0;

%% Propagation of Reference Trajectory
options = odeset('RelTol',3e-14,'AbsTol',1e-16);

[t,y] = ode45(@(t,X) jah_sat_1_ode( ...
    t, X, ENV.AccelFcn, ENV.AmatrixFcn, IC.time_struct.jd_UTC_days), ...
    IC.time.tspan, states_initial_input, options);

%% Extract States
% Extract r_sat and v_sat at the observation times from the ODE solution
X_star_table = array2table(y(:,1:7), 'VariableNames', {'r_x_km', 'r_y_km', 'r_z_km', 'v_x_km_s', 'v_y_km_s', 'v_z_km_s', 'C_drag',}); % Convert to table for easier handling
X_star_table.time_sec_past_epoch = t; % Add time column to the table

% Now add STM elements to table.
X_star_table.STM_elements = cell(height(X_star_table),1); % Pre-allocate cell array for STM elements
for i = 1:height(X_star_table)
    STM_i = reshape(y(i,8:56), N_states, N_states); % Extract STM elements for this time step
    X_star_table.STM_elements{i} = STM_i; % Store in table
end

% Now trim down to just the observation times (this will make it easier to work with for the rest of the filter)
t_obs = ENV.ref_data.Actual_Measurements.time_sec_past_epoch;
X_star_table = X_star_table(ismember(X_star_table.time_sec_past_epoch, t_obs), :); % Keep only rows where time matches observation times
% TODO: Fix this such that it doesn't rely on exact matches of time
% abs(STM_i(:,:,end) - ENV.ref_data.STM_21600_0_ref)./ENV.ref_data.STM_21600_0_ref
%% H Matrix and STM at Each Time Step
N_obs = length(t_obs);
STM_i = NaN(N_states, N_states, N_obs); % Pre-allocate 3D array for STM at each time step
H_i = NaN(2, N_states, N_obs); % Pre-allocate 3D array for H at each time step
H_mapped_epoch = NaN(2, N_states, N_obs);
residual_prefit = NaN(2, N_obs); % 2 x N_obs to make column vectors more naturally
Y_computed = NaN(2, N_obs); % 2 x N_obs to make column vectors more naturally

% Initialized Accumulated Normal Equations
Information_Matrix = zeros(N_states, N_states);
Residual_Vector = zeros(N_states, 1);

for i = 1:length(t_obs)
    
    STM_i(:,:,i) = X_star_table.STM_elements{i}; % Extract STM for this time step and store in 3D array

    % Station, Time, EOP_parameters for table row
    station_id = ENV.ref_data.Actual_Measurements.station_id(i);
    date_time_i = IC.UTC_epoch + seconds(ENV.ref_data.Actual_Measurements.time_sec_past_epoch(i));
    EOP_params_i = Tools.interpolate_EOP(date_time_i, "IERS");

    [r_station_t_i_ECI_km, v_station_t_i_ECI_km] = Tools.ECEF_ECI_Converter(IC.Stations(station_id).r_ECEF_km,zeros(3,1),date_time_i,"ECEF_to_ECI",EOP_params_i);
    H_i(:,:,i) = ENV.HmatrixFcn( ...
        X_star_table{i,1:3}(:), ... % r_sat_ECI_km
        X_star_table{i,4:6}(:), ... % v_sat_ECI_km_s 
        r_station_t_i_ECI_km(:), ...
        v_station_t_i_ECI_km(:));

    H_mapped_epoch(:,:,i) = H_i(:,:,i)*STM_i(:,:,i); % dG/dX |t_i * STM(t_i, t_0) to map from epoch to time of measurement
    
    % Compute measurement at time step i
    Y_computed(:,i) = ENV.MeasFcn( ...
        X_star_table{i,1:3}(:), ... % r_sat_ECI_km
        X_star_table{i,4:6}(:), ... % v_sat_ECI_km_s 
        r_station_t_i_ECI_km(:), ...
        v_station_t_i_ECI_km(:))';

    % Compute REsidual
    residual_prefit(:,i) = ENV.ref_data.Actual_Measurements{i, {'apparent_range_km', 'apparent_range_rate_km_s'}}(:) - Y_computed(:,i); % THESE ARE COLUMN VECTORS!

    % Note that the \eye(2) is to compute the inverse of the covariance matrix:
    R_cov_inv = IC.Stations(station_id).Covariance\eye(2); % Inverse of covariance matrix for this station
    Information_Matrix = Information_Matrix + H_mapped_epoch(:,:,i)' * R_cov_inv * H_mapped_epoch(:,:,i); % Accumulate Information Matrix
    Residual_Vector = Residual_Vector + H_mapped_epoch(:,:,i)' * R_cov_inv * residual_prefit(:,i); % Accumulate Residual Vector
end

%% Solve at Epoch
x_correction_0 = Information_Matrix\Residual_Vector; % Solve for correction at epoch
X_star_new_0 = X_star_0 + x_correction_0; % Update initial state estimate at epoch


%% Post Fit Residuals
epsilon_postfit_i = NaN(2, N_obs); % 2 x N_obs to make column vectors more naturally
for i = 1:length(t_obs)
    epsilon_postfit_i(:,i) = residual_prefit(:,i) - H_mapped_epoch(:,:,i)*x_correction_0; % Post-fit residuals at time step i
end

%% Plot Pre-fit Residuals

Measurement_Table = ENV.ref_data.Actual_Measurements;
Measurement_Table.residual_range_km = residual_prefit(1,:)';
Measurement_Table.residual_range_rate_km_s = residual_prefit(2,:)';
Measurement_Table.computed_range_km = Y_computed(1,:)';
Measurement_Table.computed_range_rate_km_s = Y_computed(2,:)';
Visuals.plot_station_residuals(Measurement_Table, {IC.Stations.Name}); hold on;

%% Now post-fit residuals
Measurement_Table.residual_range_km = epsilon_postfit_i(1,:)';
Measurement_Table.residual_range_rate_km_s = epsilon_postfit_i(2,:)';
Visuals.plot_station_residuals(Measurement_Table, {IC.Stations.Name});

%% Update ODE Run
states_initial_input = [X_star_new_0; STM_t0_t0(:)];

