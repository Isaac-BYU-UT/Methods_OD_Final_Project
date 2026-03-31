clear; clc;


%%
% Test with April 2, 2006, 00:00 UTC
UTC_test = datetime(2006, 4, 2, 0, 0, 0, 'TimeZone','UTC');
time_struct_test = Tools.ComputeTimeSystems(UTC_test);

%%

% T_UT1 = time_struct_test.t_UT1_j_centuries;
T_UT1 = (time_struct_test.jd_UTC_days - 2451545.0)/36525
% T_UT1 = 0.062491444

estimate_pos1 = Forces.Compute_Sun_Position(T_UT1)
estimate_pos2 = Forces.Vallado_sunPositionLowPrecision(time_struct_test.jd_UTC_days)
estimate_pos3 = Forces.Vallado_sunAlmanacLowPrecision(time_struct_test.jd_UTC_days)

builitin_pos = planetEphemeris(time_struct_test.jd_UTC_days, 'Earth', 'Sun')

(estimate_pos1 - builitin_pos')./builitin_pos'
(estimate_pos2 - builitin_pos')./builitin_pos'
(estimate_pos3 - builitin_pos')./builitin_pos'

%% Run Time Comparison

clear; clc;
% load('baseline_03_28_2026.mat') % Ensure Constants is available

%% Setup Test Data
UTC_test = datetime(2006, 4, 2, 0, 0, 0, 'TimeZone','UTC');
% Assuming your Tools class works as intended:
jd_UTC = 2453827.5; 
T_UT1 = (jd_UTC - 2451545.0)/36525;

%% Timing Comparison

% 1. Time your algorithm
% We pass the handle directly. timeit handles the execution.
my_alg_handle = @() Forces.Compute_Sun_Position(T_UT1);
% my_alg_handle = @() Forces.Vallado_sunPositionLowPrecision(jd_UTC);
% my_alg_handle = @() Forces.Vallado_sunAlmanacLowPrecision(jd_UTC);
t_mine = timeit(my_alg_handle);

% 2. Time the built-in function
% Note: planetEphemeris performance varies if kernels are already cached
t_builtin = timeit(@() planetEphemeris(jd_UTC, 'Earth', 'Sun'));

%% Display Results
fprintf('Analytical Algorithm: %.6f ms\n', t_mine * 1000);
fprintf('planetEphemeris:      %.6f ms\n', t_builtin * 1000);

t_builtin/t_mine
