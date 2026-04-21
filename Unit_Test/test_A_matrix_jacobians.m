function test_A_matrix_jacobians()
% Unit test to verify that the symbolic Jacobians (da/dr and da/dv) match
% numerical finite difference approximations.

% Load settings and environment
[S, ENV] = Setup.loadSettings('F', 'HW5', false, false);

% Nominal state values
r0_ECI_meters = [6990077.798814194;
                 1617465.311978378;
                 22679.810569245355];
v0_ECI_meters_per_sec = [-1675.13972506056;
                         7273.72441330686;
                         252.688512916741];
C_drag = 1.88;

% Station coordinates (not needed for A matrix, but for consistency)
stn_coords_ECEF_meters.Atoll = [-6143584, 1364250, 1033743]';
r0_ECEF_station_1_meters = stn_coords_ECEF_meters.Atoll;
v0_ECEF_station_1_meters_s = [0; 0; 0];

% Time structures
UTC.year = 2018;
UTC.month = 2;
UTC.day = 1;
UTC.hour = 5;
UTC.minute = 0;
UTC.seconds = 0.0;
UTC_date_time = datetime(UTC.year, UTC.month, UTC.day, ...
                         UTC.hour, UTC.minute, UTC.seconds, ...
                         'TimeZone', 'UTC');
time_struct = Tools.ComputeTimeSystems(UTC_date_time);

% Compute EOP, sun, moon positions
[~, ~, ~, R_ECEF_from_ECI] = Tools.ECEF_ECI_Converter(r0_ECEF_station_1_meters, v0_ECEF_station_1_meters_s, UTC_date_time, "ECEF_to_ECI", ENV.EOP_t0);
[r_sun_rel_earth_ECI_meters, ~] = Forces.Vallado_sunPositionLowPrecision(time_struct.jd_UTC_days);
[r_moon_rel_earth_ECI_meters, ~] = Forces.Vallado_moonPositionLowPrecision(time_struct.jd_UTC_days);
sat_is_illuminated = Forces.Vallado_sunLOS(r0_ECI_meters, r_sun_rel_earth_ECI_meters);

% Function handle for acceleration
accel_func = @(r, v) Forces.Compute_Total_Acceleration_ECI_m_s2(...
    r, v, C_drag, r_sun_rel_earth_ECI_meters, r_moon_rel_earth_ECI_meters, sat_is_illuminated, R_ECEF_from_ECI);

% Compute symbolic Jacobians
sym_da_dr = Forces.calc_da_dr(r0_ECI_meters, v0_ECI_meters_per_sec, C_drag, ...
                              r_sun_rel_earth_ECI_meters, r_moon_rel_earth_ECI_meters, sat_is_illuminated, R_ECEF_from_ECI);
sym_da_dv = Forces.calc_da_dv(r0_ECI_meters, v0_ECI_meters_per_sec, C_drag, ...
                              r_sun_rel_earth_ECI_meters, r_moon_rel_earth_ECI_meters, sat_is_illuminated, R_ECEF_from_ECI);

% Compute numerical Jacobians using finite differences
delta = 1e-5;  % Perturbation size
tol = 1e-10;   % Tolerance for comparison

% Numerical da/dr
num_da_dr = zeros(3, 3);
for i = 1:3
    r_plus = r0_ECI_meters;
    r_minus = r0_ECI_meters;
    r_plus(i) = r_plus(i) + delta;
    r_minus(i) = r_minus(i) - delta;
    a_plus = accel_func(r_plus, v0_ECI_meters_per_sec);
    a_minus = accel_func(r_minus, v0_ECI_meters_per_sec);
    num_da_dr(:, i) = (a_plus - a_minus) / (2 * delta);
end

% Numerical da/dv
num_da_dv = zeros(3, 3);
for i = 1:3
    v_plus = v0_ECI_meters_per_sec;
    v_minus = v0_ECI_meters_per_sec;
    v_plus(i) = v_plus(i) + delta;
    v_minus(i) = v_minus(i) - delta;
    a_plus = accel_func(r0_ECI_meters, v_plus);
    a_minus = accel_func(r0_ECI_meters, v_minus);
    num_da_dv(:, i) = (a_plus - a_minus) / (2 * delta);
end

% Compare and assert
diff_da_dr = abs(sym_da_dr - num_da_dr);
diff_da_dv = abs(sym_da_dv - num_da_dv);

max_diff_da_dr = max(diff_da_dr(:));
max_diff_da_dv = max(diff_da_dv(:));

fprintf('Max difference in da/dr: %e\n', max_diff_da_dr);
fprintf('Max difference in da/dv: %e\n', max_diff_da_dv);

% Assertions
assert(max_diff_da_dr < tol, 'Symbolic da/dr does not match numerical finite difference within tolerance.');
assert(max_diff_da_dv < tol, 'Symbolic da/dv does not match numerical finite difference within tolerance.');

fprintf('All Jacobian tests passed!\n');
end