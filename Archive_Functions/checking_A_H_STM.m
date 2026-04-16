clear; clc;
plotting_histograms = false;
%% Load Data .mat files and .dat files and EOP files for universal use
ref_data = Tools.get_reference_data();          % Includes A_t0_ref, H_tilde_t0_ref, STM_21600_0_ref, Light_Time_Truth_ref, Actual_Measurements
ECI_ECEF_Transform_Data = Tools.get_EOP_data(); % Includes table_EOP_IERS, table_EOP_Celestrak, and nut80
EOP_IERS = ECI_ECEF_Transform_Data.table_EOP_IERS;
EOP_Celestrak = ECI_ECEF_Transform_Data.table_EOP_Celestrak;

%% Load Functions
[Acceleration_Computation_func, A_matrix_computation_func] = Forces.master_dynamics();
[Measurement_Computation_func, H_matrix_computation_func] = Measurements.master_measurements();

%% EPOCH Time Structures
UTC_date_time_epoch_t0 = datetime(2018, 2, 1, 5, 0, 0.0, 'TimeZone','UTC'); % 1 Feb 2018, 05:00:00 UTC.
time_struct_t0 = Tools.ComputeTimeSystems(UTC_date_time_epoch_t0);

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

EOP_params_t0 = Tools.interpolate_EOP(UTC_date_time_epoch_t0, "IERS"); % Interpolate EOP for initial t_0

Stations = struct();

for i = 1:length(station_names)
    Stations(i).Name = station_names{i};
    Stations(i).r_ECEF_km = raw_coords_ECEF_m(i, :)' * Constants.METERS_TO_KM;
    
    [r_t0_ECI_km, v_t0_ECI_km_s] = Tools.ECEF_ECI_Converter(Stations(i).r_ECEF_km, zeros(3,1), ...
                                                                UTC_date_time_epoch_t0, 'ECEF_to_ECI', EOP_params_t0);
    
    Stations(i).r_t0_ECI_km =  r_t0_ECI_km;
    Stations(i).v_t0_ECI_km_s = v_t0_ECI_km_s;
end

%% Compute Initial Sun and Moon Positions in ECI at t_0
% [r_sun_t0_rel_earth_ECI_km, ~] = planetEphemeris(time_struct_t0.jd_UTC_days, 'Earth', 'Sun');
% [r_moon_t0_rel_earth_ECI_km, ~] = planetEphemeris(time_struct_t0.jd_UTC_days, 'Earth', 'Moon');

[r_sun_t0_rel_earth_ECI_km, ~] = Forces.Vallado_sunPositionLowPrecision(time_struct_t0.jd_UTC_days);
[r_moon_t0_rel_earth_ECI_km, ~] = Forces.Vallado_moonPositionLowPrecision(time_struct_t0.jd_UTC_days);
sat_is_illuminated = Forces.Vallado_sunLOS(r_sat_t0_ECI_km,r_sun_t0_rel_earth_ECI_km);

%% Question 1 (A, H matrix)

A_t0 = A_matrix_computation_func(r_sat_t0_ECI_km, v_sat_t0_ECI_km_s, C_drag, r_sun_t0_rel_earth_ECI_km, r_moon_t0_rel_earth_ECI_km, sat_is_illuminated);
H_tilde_t0 = H_matrix_computation_func(r_sat_t0_ECI_km, v_sat_t0_ECI_km_s, Stations(1).r_t0_ECI_km, Stations(1).v_t0_ECI_km_s); % Use station 1 for this comparison

%% Question 2 (STM)

% Time Vector
time_ode.end = 21600; % sec 
time_ode.store_interval = 60; %sec
time_ode.tspan = 0:time_ode.store_interval:time_ode.end; % sec

% Initial Conditions
STM_0 = eye(7); y0 = [r_sat_t0_ECI_km; v_sat_t0_ECI_km_s; C_drag; STM_0(:)];

% Run Integration
options = odeset('RelTol', 3e-14, 'AbsTol', 1e-16);
[t, y] = ode45(@(t,X) jah_sat_1_ode(t, X, Acceleration_Computation_func, A_matrix_computation_func, time_struct_t0.jd_UTC_days),time_ode.tspan, y0, options);

%% Compare to Reference Data

relDiff_A = abs((A_t0 - ref_data.A_t0_ref)./ref_data.A_t0_ref) % Improve with better pertubation models
relDiff_H = abs((H_tilde_t0 - ref_data.H_tilde_t0_ref)./ref_data.H_tilde_t0_ref) % Improve with light-time corrections
STM_21600_0 = reshape(y(end,8:56),7,7);
relDiff_STM = abs((STM_21600_0 - ref_data.STM_21600_0_ref)./ref_data.STM_21600_0_ref)

if(plotting_histograms) 
    Visuals.plot_error_exponents(relDiff_A);
    Visuals.plot_error_exponents(relDiff_H);
    Visuals.plot_error_exponents(relDiff_STM); 
end

%% Question 3

% Extract r_sat and v_sat at the observation times from the ODE solution
table_satellite_states = array2table(y(:,1:7), 'VariableNames', {'r_x_km', 'r_y_km', 'r_z_km', 'v_x_km_s', 'v_y_km_s', 'v_z_km_s', 'C_drag'}); % Convert to table for easier handling
table_satellite_states.time_sec_past_epoch = t; % Add time column to the table

Measurement_Comp_Table = ref_data.Actual_Measurements;
Measurement_Comp_Table.computed_range_km = NaN(height(Measurement_Comp_Table),1);        % Add this column to populate
Measurement_Comp_Table.computed_range_rate_km_s = NaN(height(Measurement_Comp_Table),1); % Add this column to populate

for i = 1:height(Measurement_Comp_Table)

    % Station, Time, EOP_parameters for table row
    station_id = Measurement_Comp_Table.station_id(i);
    date_time_i = UTC_date_time_epoch_t0 + seconds(Measurement_Comp_Table.time_sec_past_epoch(i));
    EOP_params_i = Tools.interpolate_EOP(date_time_i, "IERS");

    % Find Computed State Vector at Time
    j = find(t == Measurement_Comp_Table.time_sec_past_epoch(i),1);
    j_state_output = y(j, :);
    r_sat_ti_ECI_km = j_state_output(1:3);
    v_sat_ti_ECI_km_s = j_state_output(4:6);

    % Compute Station Coordinates at State, Range/Range-Rate, & Store
    [r_station_t_i_ECI_km, v_station_t_i_ECI_km] = Tools.ECEF_ECI_Converter(Stations(station_id).r_ECEF_km,zeros(3,1),date_time_i,"ECEF_to_ECI",EOP_params_i);
    range_range_rate_output = Measurement_Computation_func(r_sat_ti_ECI_km(:),v_sat_ti_ECI_km_s(:),r_station_t_i_ECI_km, v_station_t_i_ECI_km);
    Measurement_Comp_Table.computed_range_km(i) = range_range_rate_output(1);
    Measurement_Comp_Table.computed_range_rate_km_s(i) = range_range_rate_output(2);
end

% ... After your loop that populates Measurement_Comp_Table ...

% 1. Create a list of names for the legend
station_list = {Stations.Name}; 

% 2. Run Visualization Functions
Visuals.plot_station_residuals(Measurement_Comp_Table, station_list);
Visuals.plot_measurement_correlation(Measurement_Comp_Table);
Visuals.plot_measurement_correlation_linked(Measurement_Comp_Table);

% 3. Bonus: Simple Histogram to check for Zero-Mean Residuals (Gausianity)
range_residuals_km = Measurement_Comp_Table.apparent_range_km - Measurement_Comp_Table.computed_range_km;
range_rate_residuals_km = Measurement_Comp_Table.apparent_range_rate_km_s - Measurement_Comp_Table.computed_range_rate_km_s;

figure('Color', 'w');
histogram(range_residuals_km, 20, 'FaceColor', [0.4 0.4 0.4]);
xlabel('Residual Error (km)'); ylabel('Frequency');
title('Range Residual Distribution (Should be Zero-Centered)');

%% Plot Orbit
r_final = y(end,1:3);
v_final = y(end,4:6);
full_orbit_ECI = y(:,1:6); % Extract full orbit states for plotting
Visuals.plot_orbit_errors(r_final(:),v_final(:),eye(6)) % TODO: Figure out weird covariance stuff
Visuals.plot_position(r_final(:),v_final(:),eye(6), full_orbit_ECI)


%% Figuring Out Light Time Stuff
figure;
yyaxis left;
plot(ref_data.Light_Time_Truth_ref.time_since_epoch_sec,ref_data.Light_Time_Truth_ref.r_z_km,'bx')
yyaxis right;
plot(t,y(:,3),'rx')
