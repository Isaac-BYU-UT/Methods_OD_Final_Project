clear; clc;

%% Load Initial Conditions
IC = Initial_Conditions.initial_conditions_HW5();

%% Setup Environment
ENV = Environment.setup_environment(IC);

%%

%  ODE INITIAL VECTOR
STM_t0_t0 = eye(7);
IC.y0 = [IC.r_sat_ECI_km;
         IC.v_sat_ECI_km_s;
         IC.C_drag;
         STM_t0_t0(:)];

%% A and H
A_t0 = ENV.AmatrixFcn( ...
    IC.r_sat_ECI_km, ...
    IC.v_sat_ECI_km_s, ...
    IC.C_drag, ...
    ENV.r_sun_ECI_km, ...
    ENV.r_moon_ECI_km, ...
    ENV.sunLOS);

H_tilde_t0 = ENV.HmatrixFcn( ...
                IC.r_sat_ECI_km, ...
                IC.v_sat_ECI_km_s, ...
                ENV.Stations(1).r_ECI_km, ...
                ENV.Stations(1).v_ECI_km_s);

%% Propagation
options = odeset('RelTol',3e-14,'AbsTol',1e-16);

[t,y] = ode45(@(t,X) jah_sat_1_ode( ...
    t, X, ENV.AccelFcn, ENV.AmatrixFcn, IC.time_struct.jd_UTC_days), ...
    IC.time.tspan, IC.y0, options);

%% Compare to Reference Data

relDiff_A = abs((A_t0 - ENV.ref_data.A_t0_ref)./ENV.ref_data.A_t0_ref) % Improve with better pertubation models
relDiff_H = abs((H_tilde_t0 - ENV.ref_data.H_tilde_t0_ref)./ENV.ref_data.H_tilde_t0_ref) % Improve with light-time corrections
STM_21600_0 = reshape(y(end,8:56),7,7);
relDiff_STM = abs((STM_21600_0 - ENV.ref_data.STM_21600_0_ref)./ENV.ref_data.STM_21600_0_ref)

if(IC.plotting_histograms) 
    Visuals.plot_error_exponents(relDiff_A);
    Visuals.plot_error_exponents(relDiff_H);
    Visuals.plot_error_exponents(relDiff_STM); 
end

%% Question 3

% Extract r_sat and v_sat at the observation times from the ODE solution
table_satellite_states = array2table(y(:,1:7), 'VariableNames', {'r_x_km', 'r_y_km', 'r_z_km', 'v_x_km_s', 'v_y_km_s', 'v_z_km_s', 'C_drag'}); % Convert to table for easier handling
table_satellite_states.time_sec_past_epoch = t; % Add time column to the table

Measurement_Comp_Table = ENV.ref_data.Actual_Measurements;
Measurement_Comp_Table.computed_range_km = NaN(height(Measurement_Comp_Table),1);        % Add this column to populate
Measurement_Comp_Table.computed_range_rate_km_s = NaN(height(Measurement_Comp_Table),1); % Add this column to populate

for i = 1:height(Measurement_Comp_Table)

    % Station, Time, EOP_parameters for table row
    station_id = Measurement_Comp_Table.station_id(i);
    date_time_i = IC.UTC_epoch + seconds(Measurement_Comp_Table.time_sec_past_epoch(i));
    EOP_params_i = Tools.interpolate_EOP(date_time_i, "IERS");

    % Find Computed State Vector at Time
    j = find(t == Measurement_Comp_Table.time_sec_past_epoch(i),1);
    j_state_output = y(j, :);
    r_sat_ti_ECI_km = j_state_output(1:3);
    v_sat_ti_ECI_km_s = j_state_output(4:6);

    % Compute Station Coordinates at State, Range/Range-Rate, & Store
    [r_station_t_i_ECI_km, v_station_t_i_ECI_km] = Tools.ECEF_ECI_Converter(IC.Stations(station_id).r_ECEF_km,zeros(3,1),date_time_i,"ECEF_to_ECI",EOP_params_i);
    range_range_rate_output = ENV.MeasFcn(r_sat_ti_ECI_km(:),v_sat_ti_ECI_km_s(:),r_station_t_i_ECI_km, v_station_t_i_ECI_km);
    Measurement_Comp_Table.computed_range_km(i) = range_range_rate_output(1);
    Measurement_Comp_Table.computed_range_rate_km_s(i) = range_range_rate_output(2);
end

range_residuals_km = Measurement_Comp_Table.apparent_range_km - Measurement_Comp_Table.computed_range_km;
range_rate_residuals_km = Measurement_Comp_Table.apparent_range_rate_km_s - Measurement_Comp_Table.computed_range_rate_km_s;


% Calculate RMS and Display
range_residual_RMS_km = sqrt(mean(range_residuals_km.^2));
range_rate_residual_RMS_km_s = sqrt(mean(range_rate_residuals_km.^2));
disp(['Range Residual RMS: ', num2str(range_residual_RMS_km), ' km']);
disp(['Range-rate Residual RMS: ', num2str(range_rate_residual_RMS_km_s), ' km/s']);

% 1. Create a list of names for the legend
station_list = {IC.Stations.Name}; 

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
plot(ENV.ref_data.Light_Time_Truth_ref.time_since_epoch_sec,ENV.ref_data.Light_Time_Truth_ref.r_z_km,'bx')
yyaxis right;
plot(t,y(:,3),'rx')