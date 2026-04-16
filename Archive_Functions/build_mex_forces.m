% build_mex_forces.m

% 1. Define typical input types (Crucial for C static typing)
% All inputs must be doubles to match symbolic output
r_type   = coder.typeof(double(0), [3 1], false); % 3x1 fixed-size double
v_type   = coder.typeof(double(0), [3 1], false); % 3x1 fixed-size double
Cd_type  = double(0);                             % Scalar double
sun_type = coder.typeof(double(0), [3 1], false);
mon_type = coder.typeof(double(0), [3 1], false);
ill_type = double(0);                             % Logical/Scalar double
R_type   = coder.typeof(double(0), [3 3], false); % 3x3 fixed-size double

% Combine into a cell array for the -args flag
input_args = {r_type, v_type, Cd_type, sun_type, mon_type, ill_type, R_type};

% 2. Setup Coder Configuration for Speed
cfg = coder.config('mex');
cfg.IntegrityChecks = false; % Big speed boost: disables array bounds checking
cfg.ResponsivenessChecks = false; % Big speed boost: allows the MEX to run without interruption

% 3. Run Codegen
% Note: We call 'Forces.function' but output to a plain name to avoid '+' issues
fprintf('Compiling MEX files...\n');

codegen -config cfg Forces.calc_da_dr  -args input_args -o calc_da_dr_mex
codegen -config cfg Forces.calc_da_dv  -args input_args -o calc_da_dv_mex
codegen -config cfg Forces.calc_da_dCd -args input_args -o calc_da_dCd_mex
codegen -config cfg Forces.Compute_Total_Acceleration_ECI_m_s2 -args input_args -o accel_mex

fprintf('Done! Use calc_da_dr_mex(...) in your loop.\n');