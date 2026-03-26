clear; clc;
%% Load Data
% Load .mat files and .dat files and EOP files for universal use
ref_data = get_reference_data(); % Includes A_t0_ref, H_tilde_t0_ref, STM_21600_0_ref, Light_Time_Truth_ref, Actual_Measurements
ECI_ECEF_Transform_Data = get_EOP_data(); % Includes table_EOP_IERS, table_EOP_Celestrak, and nut80
EOP_IERS = ECI_ECEF_Transform_Data.table_EOP_IERS;
EOP_Celestrak = ECI_ECEF_Transform_Data.table_EOP_Celestrak;

%% Load Functions
[Acceleration_Computation_func, A_matrix_computation_func] = Forces.master_dynamics();
[Measurement_Computation_func, H_matrix_computation_func] = Measurements.master_measurements();

%% Time Structures

% Gregorian Date (UTC) 1 Feb 2018, 05:00:00 UTC.
UTC.year = 2018; UTC.month = 2; UTC.day = 1; UTC.hour = 5; UTC.minute = 0; UTC.seconds = 0.0;
UTC_date_time_epoch_t0 = datetime(UTC.year, UTC.month, UTC.day,UTC.hour, UTC.minute, UTC.seconds, 'TimeZone','UTC');

%% Satellite Initial Coordinates
r_sat_t0_ECI_km   = [6990077.798814194;1617465.311978378;22679.810569245355] * Constants.METERS_TO_KM;
v_sat_t0_ECI_km_s = [-1675.13972506056;7273.72441330686;252.688512916741] * Constants.METERS_TO_KM;
C_drag = 1.88;

%% Station Initial Coordinates --> ECI

% Define names and coordinates in a clean list
station_names = {'Atoll', 'Diego Garcia', 'Arecibo'};
raw_coords_ECEF_m = [
    -6143584,  1364250,  1033743;  % Atoll
     1907295,  6030810,  -817119;  % Diego Garcia
     2390310, -5564341,  1994578   % Arecibo
];

% Interpolate EOP once
EOP_params_t0 = interpolate_EOP(UTC_date_time_epoch_t0, "IERS");

% Initialize the structure
Stations = struct();

for i = 1:length(station_names)
    Stations(i).Name = station_names{i};
    Stations(i).r_ECEF_km = raw_coords_ECEF_m(i, :)' * Constants.METERS_TO_KM;
    
    % Convert to ECI in one loop
    [r_t0_ECI_km, v_t0_ECI_km_s] = ECEF_ECI_Converter(Stations(i).r_ECEF_km, zeros(3,1), ...
                        UTC_date_time_epoch_t0, 'ECEF_to_ECI', EOP_params_t0);
    
    Stations(i).r_t0_ECI_km =  r_t0_ECI_km;
    Stations(i).v_t0_ECI_km_s = v_t0_ECI_km_s;
end

%% Question 1
% Derive the An× n and ˜Hm× n 
% matrices for the linearized system and implement the partials in your
% computer language of choice. We recommend using a symbolic solver, e.g., MATLAB’s symbolic
% toolbox, due to the large number of partial derivatives required. Compare to the numeric solutions
% online (found on canvas) at t0. Compute the relative difference for each non-zero element.
%
% For the ˜H solution, provide the numeric values for the relative difference in your write-up. For
% the A matrix, include a histogram of the exponents.
A_t0 = A_matrix_computation_func(r_sat_t0_ECI_km, v_sat_t0_ECI_km_s, C_drag);
H_tilde_t0 = H_matrix_computation_func(r_sat_t0_ECI_km, v_sat_t0_ECI_km_s, Stations(1).r_t0_ECI_km, Stations(1).v_t0_ECI_km_s); % Use station 1 for this comparison

relDiff_A = abs((A_t0 - ref_data.A_t0_ref)./ref_data.A_t0_ref) % Improve with better pertubation models
relDiff_H = abs((H_tilde_t0 - ref_data.H_tilde_t0_ref)./ref_data.H_tilde_t0_ref) % Improve with light-time corrections

plot_error_exponents(relDiff_A);
plot_error_exponents(relDiff_H);

%% Question 2
% Compare STM at time t = 21600 seconds
% Integrate position, velocity, and Φ(ti,t0) 
% from t = 0,...,21600 seconds. Store the results in 60
% second intervals. Compare your results with those on the web (found on the Canvas site) and
% compute the relative difference of the top left relevant portion of Φ(21600,0). Like the previous
% question, provide a histogram of the exponents for the STM comparison.

% Time Vector
time_ode.start = 0; % sec
time_ode.end = 21600; % sec
time_ode.integration_interval = .1; % sec
time_ode.store_interval = 60; %sec
time_ode.tspan = time_ode.start:time_ode.integration_interval:time_ode.end; % sec

% Initial Conditions
STM_0 = eye(7); y0 = [r_sat_t0_ECI_km; v_sat_t0_ECI_km_s; C_drag; STM_0(:)];

% Run Integration
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-16);
[t, y] = ode45(@(t,X) jah_sat_1_ode(t, X, Acceleration_Computation_func, A_matrix_computation_func),time_ode.tspan, y0, options);

% Extract data from t and y at store_interval
t_store = t(mod(t, time_ode.store_interval) == 0);
y_store = y(mod(t, time_ode.store_interval) == 0, :);

STM_21600_0 = reshape(y(end,8:56),7,7);
relDiff_STM = abs((STM_21600_0 - ref_data.STM_21600_0_ref)./ref_data.STM_21600_0_ref)
plot_error_exponents(relDiff_STM);

%% Question 3
% Calculate the predicted range and range-rate for the appropriate tracking station at each observa
% tion time (use the data provided on the Canvas site). Compare and plot the range and range-rate
% residuals (post-fit). Calculate the range residual RMS and range-rate residual RMS

% Extract r_sat and v_sat at the observation times from the ODE solution
table_satellite_states = array2table(y(:,1:7), 'VariableNames', {'r_x_km', 'r_y_km', 'r_z_km', 'v_x_km_s', 'v_y_km_s', 'v_z_km_s', 'C_drag'}); % Convert to table for easier handling
table_satellite_states.time_sec_past_epoch = t; % Add time column to the table

% Now trim down the satellite states to only the observation times (every 60 seconds)
table_satellite_states_observation_times = table_satellite_states(mod(table_satellite_states.time_sec_past_epoch, time_ode.store_interval) == 0, :);
table_satellite_states_observation_times.date_time_UTC = UTC_date_time_epoch_t0 + seconds(table_satellite_states_observation_times.time_sec_past_epoch);

Measurement_Comp_Table = ref_data.Actual_Measurements;
Measurement_Comp_Table.computed_range_km = NaN(height(Measurement_Comp_Table),1);
Measurement_Comp_Table.computed_range_rate_km_s = NaN(height(Measurement_Comp_Table),1);

for i = 1:height(Measurement_Comp_Table)
    % Compute EOP for time step:
    station_id = Measurement_Comp_Table.station_id(i);
    date_time_i = UTC_date_time_epoch_t0 + seconds(Measurement_Comp_Table.time_sec_past_epoch(i));
    EOP_params_i = interpolate_EOP(date_time_i, "IERS");
    [r_station_t_i_ECI_km, v_station_t_i_ECI_km] = ECEF_ECI_Converter(Stations(station_id).r_ECEF_km,zeros(3,1),date_time_i,"ECEF_to_ECI",EOP_params_i);
    j = find(t_store == Measurement_Comp_Table.time_sec_past_epoch(i),1);
    j_state_output = y_store(j, :);
    r_sat_ti_ECI_km = j_state_output(1:3);
    v_sat_ti_ECI_km_s = j_state_output(4:6);
    range_range_rate_output = Measurement_Computation_func(r_sat_ti_ECI_km(:),v_sat_ti_ECI_km_s(:),r_station_t_i_ECI_km, v_station_t_i_ECI_km);
    Measurement_Comp_Table.computed_range_km(i) = range_range_rate_output(1);
    Measurement_Comp_Table.computed_range_rate_km_s(i) = range_range_rate_output(2);
end

range_residuals_km = Measurement_Comp_Table.apparent_range_km - Measurement_Comp_Table.computed_range_km;
range_rate_residuals_km = Measurement_Comp_Table.apparent_range_rate_km_s - Measurement_Comp_Table.computed_range_rate_km_s;

% First figure: range residuals
figure;
% Added 'MarkerSize' and 'LineWidth'
plot(Measurement_Comp_Table.time_sec_past_epoch, range_residuals_km, 'x', ...
     'MarkerSize', 8, 'LineWidth', 1.2); 
xlabel('Time since epoch (s)');
ylabel('Range residual (km)');
title('Range Residuals vs Time');
grid on;
set(gca,'FontSize', 16);

% Second figure: range-rate residuals
figure;
% Added 'MarkerSize' and 'LineWidth'
plot(Measurement_Comp_Table.time_sec_past_epoch, range_rate_residuals_km, 'x', ...
     'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('Time since epoch (s)');
ylabel('Range-rate residual (km/s)');
title('Range-rate Residuals vs Time');
grid on;
set(gca,'FontSize', 16);

%% Helper Functions --> To move to a separate file.


