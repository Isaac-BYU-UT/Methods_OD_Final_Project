clear; clc;

%% ===================== SETUP =====================
[S, ENV] = Setup.loadSettings('F','1 Hour','Final - 1 Day Subset',false,true);

N_states = 7;

% Initial STM and State
STM_0 = eye(N_states);
X_nominal_0 = [ ...
    S.IC_Sat_Epoch.position_ECI_meters;
    S.IC_Sat_Epoch.velocity_ECI_meters_per_second;
    1.88;
    STM_0(:)
];

% Time setup
time_struct = Tools.ComputeTimeSystems(S.IC_Sat_Epoch.epoch_date_time_UTC);
ENV.EOP_IERS.MJD_days = time_struct.mjd_UTC_days;

% Integrator options
options = odeset('RelTol',3e-10,'AbsTol',1e-12);

% A priori covariance
P_cov = S.StateCovariances;

% Observations
t_obs = S.ref_data.Actual_Measurements.time_sec_past_epoch;
N_obs = length(t_obs);

% Storage
Y_prefit  = zeros(2, N_obs);
Y_postfit = zeros(2, N_obs);

%% ===================== INITIALIZATION =====================
X_input = X_nominal_0;
t_prev  = t_obs(1);

%% ===================== FIRST OBSERVATION =====================
if t_obs(1) == 0

    station_id = S.ref_data.Actual_Measurements.station_id(1);

    y_meas = S.ref_data.Actual_Measurements{1, ...
        {'apparent_range_meters','apparent_range_rate_meters_s'}}(:);

    R = S.Stations(station_id).Covariance;

    % Station state
    date_time = S.IC_Sat_Epoch.epoch_date_time_UTC;
    EOP = Tools.interpolate_EOP(date_time, ENV.EOP_IERS, ENV.EOP_Celestrak);

    [r_stn, v_stn] = Tools.ECEF_ECI_Converter( ...
        S.Stations(station_id).position_ECEF_meters, ...
        zeros(3,1), date_time, "ECEF_to_ECI", EOP);

    % ---- PREFIT ----
    X_nominal = X_input(1:N_states);

    y_computed = Measurements.Compute_Range_Range_Rate( ...
        X_nominal(1:3), X_nominal(4:6), r_stn, v_stn);

    Y_prefit(:,1) = y_computed;
    residual = y_meas - y_computed;

    % ---- UPDATE CHECK ----
    update = strcmp(active_station_mode,'All') || ...
             station_id == station_map.(active_station_mode);

    if update
        H = ENV.HmatrixFcn(X_nominal(1:3), X_nominal(4:6), r_stn, v_stn);

        K = P_cov * H' / (H * P_cov * H' + R);

        dx = K * residual;
        X_nominal = X_nominal + dx;

        I = eye(N_states);
        P_cov = (I - K*H)*P_cov*(I - K*H)' + K*R*K';
    end

    % ---- POSTFIT ----
    y_post = Measurements.Compute_Range_Range_Rate( ...
        X_nominal(1:3), X_nominal(4:6), r_stn, v_stn);

    Y_postfit(:,1) = y_post;

    % Reset STM
    X_input = [X_nominal; STM_0(:)];

    fprintf('t=0 processed. Update applied: %d\n', update);
end

%% ===================== MAIN EKF LOOP =====================
for i = 2:N_obs

    %% ---- LOAD OBSERVATION ----
    station_id = S.ref_data.Actual_Measurements.station_id(i);

    y_meas = S.ref_data.Actual_Measurements{i, ...
        {'apparent_range_meters','apparent_range_rate_meters_s'}}(:);

    R = S.Stations(station_id).Covariance;

    %% ---- PROPAGATE STATE ----
    [t, y] = ode45(@(t,X) jah_sat_1_ode(t, X, ...
        time_struct.jd_UTC_days, false), ...
        [t_obs(i-1), t_obs(i)], X_input, options);

    X_full = y(end,:)';
    X_nominal = X_full(1:N_states);

    STM = reshape(X_full(N_states+1:end), N_states, N_states);

    r_sat = X_nominal(1:3);
    v_sat = X_nominal(4:6);

    %% ---- TIME UPDATE (COVARIANCE) ----
    Q = zeros(N_states); % Process noise OFF for now

    P_bar = STM * P_cov * STM' + Q;

    %% ---- STATION STATE ----
    date_time = time_struct.jd_UTC_days + seconds(t_obs(i));

    EOP = Tools.interpolate_EOP(date_time, ENV.EOP_IERS, ENV.EOP_Celestrak);

    [r_stn, v_stn] = Tools.ECEF_ECI_Converter( ...
        S.Stations(station_id).position_ECEF_meters, ...
        zeros(3,1), date_time, "ECEF_to_ECI", EOP);

    %% ---- PREFIT ----
    y_computed = Measurements.Compute_Range_Range_Rate( ...
        r_sat, v_sat, r_stn, v_stn);

    Y_prefit(:,i) = y_computed;
    residual = y_meas - y_computed;

    %% ---- UPDATE CHECK ----
    update = strcmp(active_station_mode,'All') || ...
             station_id == station_map.(active_station_mode);

    if update
        H = ENV.HmatrixFcn(r_sat, v_sat, r_stn, v_stn);

        K = P_bar * H' / (H * P_bar * H' + R);

        dx = K * residual;
        X_nominal = X_nominal + dx;

        I = eye(N_states);
        P_cov = (I - K*H)*P_bar*(I - K*H)' + K*R*K';
        P_cov = 0.5 * (P_cov + P_cov'); % enforce symmetry
    else
        P_cov = P_bar;
        dx = zeros(N_states,1);
    end

    %% ---- POSTFIT ----
    y_post = Measurements.Compute_Range_Range_Rate( ...
        X_nominal(1:3), X_nominal(4:6), r_stn, v_stn);

    Y_postfit(:,i) = y_post;

    %% ---- LOGGING ----
    if print_updates || mod(i,20)==0
        fprintf('Step %d / %d\n', i, N_obs);
        disp('State correction:'); disp(dx');
        disp('Covariance:'); disp(P_cov);
    end

    %% ---- PREP NEXT STEP ----
    X_input = [X_nominal; eye(N_states)(:)];
end