clear; clc;
[S,ENV] = Setup.loadSettings('F','HW5',true,true);
%% Satellite Initial Coordinates
r0_ECI_meters =     [6990077.798814194;
                    1617465.311978378; 
                    22679.810569245355];
v0_ECI_meters_per_sec =     [-1675.13972506056; 
                             7273.72441330686; 
                             252.688512916741];

R0_ECI_meters = r0_ECI_meters;
V0_ECI_meters_s = v0_ECI_meters_per_sec;
C_drag = 1.88;

%% Station Initial Coordinates --> ECI

% 1 Feb 2018, 05:00:00 UTC
% Build Time Structures

%  Initial Satellite States:
stn_coords_ECEF_meters.Atoll       = [-6143584, 1364250, 1033743]'; % Make sure these are column vectors!
stn_coords_ECEF_meters.DiegoGarcia = [1907295, 6030810, -817119]';
stn_coords_ECEF_meters.Arecibo     = [2390310, -5564341, 1994578]';

r0_ECEF_station_1_meters = stn_coords_ECEF_meters.Atoll;
v0_ECEF_station_1_meters_s = [0;0;0]; % Station not moving relative to ECEF

r0_ECEF_station_2_meters = stn_coords_ECEF_meters.DiegoGarcia;
v0_ECEF_station_2_meters_s = [0;0;0]; % Station not moving relative to ECEF

r0_ECEF_station_3_meters = stn_coords_ECEF_meters.Arecibo;
v0_ECEF_station_3_meters_s = [0;0;0]; % Station not moving relative to ECEF

%% Time Structures

% Gregorian Date (UTC) 1 Feb 2018, 05:00:00 UTC.
UTC.year = 2018;
UTC.month = 2;
UTC.day = 1;
UTC.hour = 5;
UTC.minute = 0;
UTC.seconds = 0.0;

UTC_date_time = datetime(UTC.year, UTC.month, UTC.day, ...
                   UTC.hour, UTC.minute, UTC.seconds, ...
                   'TimeZone','UTC');

time_struct = Tools.ComputeTimeSystems(UTC_date_time);
%% Find EOP Params
[r0_ECI_station1_meters, v0_ECI_station1_meters_s,~,R_ECEF_from_ECI] = Tools.ECEF_ECI_Converter(r0_ECEF_station_1_meters, v0_ECEF_station_1_meters_s, UTC_date_time, "ECEF_to_ECI", ENV.EOP_t0);

[r_sun_rel_earth_ECI_meters, ~] = Forces.Vallado_sunPositionLowPrecision(time_struct.jd_UTC_days);
[r_moon_rel_earth_ECI_meters, ~] = Forces.Vallado_moonPositionLowPrecision(time_struct.jd_UTC_days);
sat_is_illuminated = Forces.Vallado_sunLOS(r0_ECI_meters,r_sun_rel_earth_ECI_meters);
%% Compute A and H Matrices

A_t0 = Forces.get_A_matrix(R0_ECI_meters, V0_ECI_meters_s,C_drag,r_sun_rel_earth_ECI_meters,r_moon_rel_earth_ECI_meters,sat_is_illuminated,R_ECEF_from_ECI);
H_tilde_t0 = Measurements.Compute_H_matrix(R0_ECI_meters, V0_ECI_meters_s, r0_ECI_station1_meters, v0_ECI_station1_meters_s); % Station 1 is the measurement for t = 0.

A_t0_ref = S.ref_data.A_t0_ref;
H_tilde_t0_ref = S.ref_data.H_tilde_t0_ref;

relDiff_A = abs((A_t0 - A_t0_ref)./A_t0_ref)
relDiff_H = abs((H_tilde_t0 - H_tilde_t0_ref)./H_tilde_t0_ref)

%% Propogate State/STM

% Time Vector
t_start = 0; % sec
t_end = 21600; % sec
t_interval_ode = 1; % sec
t_interval_store = 60; %sec
tspan = t_start: t_interval_ode : t_end;

% Options
relative_tolerance = 3e-14;           % Same tolerances dictated in HW02
abs_tolerance = 1e-16;
options = odeset('RelTol', relative_tolerance, 'AbsTol', abs_tolerance);

% Initial Conditions
STM_0 = eye(7);
y0 = [R0_ECI_meters; V0_ECI_meters_s; C_drag; STM_0(:)];

% Run Integration
[t, y] = ode45(@(t,X) jah_sat_1_ode(t, X, S, ENV, false), ...
               tspan, y0, options);


STM_21600_0_ref = S.ref_data.STM_21600_0_ref;

STM_21600_0 = reshape(y(end,8:56),7,7);

abs((STM_21600_0 - STM_21600_0_ref)./STM_21600_0_ref)






