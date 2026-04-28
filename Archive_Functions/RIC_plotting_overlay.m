% RSW Covariance Projection and Plotting Script (REVISED + SW)
clear; clc; close all;

%% 1. Define Files and Reference
results_dir = 'Results/';
ref_filename = 'EKF_Results_All_station_All_data_2026-04-02_11-17-14.mat';

files_to_plot = {
    'EKF_Results_All_station_All_data_2026-04-02_11-17-14.mat',       'All Stations, All Data',      'k';
    'EKF_Results_Atoll_station_All_data_2026-04-02_11-20-19.mat',     'Atoll Station, All Data',     'r';
    'EKF_Results_Diego_Garcia_station_All_data_2026-04-02_11-23-18.mat','Diego Garcia, All Data',  'b';
    'EKF_Results_All_station_Range_data_2026-04-02_11-36-28.mat',     'All Stations, Range Only',    'm';
    'EKF_Results_All_station_Range_Rate_data_2026-04-02_11-40-38.mat','All Stations, Rate Only',     'c'
};

%% 2. Load Reference State
ref_path = fullfile(results_dir, ref_filename);

ref_data = load(ref_path, ...
    'r_sat_t_i_ECI_km', ...
    'v_sat_t_i_ECI_km_s');

r_ref_ECI = ref_data.r_sat_t_i_ECI_km;
v_ref_ECI = ref_data.v_sat_t_i_ECI_km_s;

%% 3. Prepare Figures
fig_RI = figure('Name', 'Radial vs Intrack', 'Position', [100,100,800,600]);
hold on; grid on; axis equal;
xlabel('Intrack Error (km)', 'FontWeight', 'bold');
ylabel('Radial Error (km)', 'FontWeight', 'bold');
title('3\sigma Position Uncertainty (Radial vs Intrack)');

fig_RC = figure('Name', 'Radial vs Crosstrack', 'Position', [150,150,800,600]);
hold on; grid on; axis equal;
xlabel('Crosstrack Error (km)', 'FontWeight', 'bold');
ylabel('Radial Error (km)', 'FontWeight', 'bold');
title('3\sigma Position Uncertainty (Radial vs Crosstrack)');

% ✅ NEW FIGURE
fig_SW = figure('Name', 'Intrack vs Crosstrack', 'Position', [200,200,800,600]);
hold on; grid on; axis equal;
xlabel('Intrack Error (km)', 'FontWeight', 'bold');
ylabel('Crosstrack Error (km)', 'FontWeight', 'bold');
title('3\sigma Position Uncertainty (Intrack vs Crosstrack)');

handles_RI = []; labels_RI = {};
handles_RC = []; labels_RC = {};
handles_SW = []; labels_SW = {};

%% 4. Loop Through Cases
for i = 1:size(files_to_plot,1)

    filename = files_to_plot{i,1};
    label    = files_to_plot{i,2};
    color    = files_to_plot{i,3};

    filepath = fullfile(results_dir, filename);

    try
        data = load(filepath, ...
            'r_sat_t_i_ECI_km', ...
            'v_sat_t_i_ECI_km_s', ...
            'P_cov_i_minus_1');
    catch
        warning('Skipping %s (missing data)', filename);
        continue;
    end

    r_case_ECI = data.r_sat_t_i_ECI_km;
    v_case_ECI = data.v_sat_t_i_ECI_km_s;

    P_full = data.P_cov_i_minus_1;

    if ndims(P_full) == 3
        P_ECI = P_full(1:3,1:3,end);
    else
        P_ECI = P_full(1:3,1:3);
    end

    % --- RSW Frame ---
    R_unit = r_case_ECI / norm(r_case_ECI);
    W_unit = cross(r_case_ECI, v_case_ECI); W_unit = W_unit / norm(W_unit);
    S_unit = cross(W_unit, R_unit);

    T = [R_unit'; S_unit'; W_unit'];

    % --- Transform ---
    delta_r_ECI = r_case_ECI - r_ref_ECI;
    delta_r_RSW = T * delta_r_ECI;

    P_RSW = T * P_ECI * T';
    P_RSW = (P_RSW + P_RSW') / 2;

    R_err = delta_r_RSW(1);
    S_err = delta_r_RSW(2);
    W_err = delta_r_RSW(3);

    % --- 2D Covariances ---
    P_RI = [P_RSW(2,2), P_RSW(2,1);
            P_RSW(1,2), P_RSW(1,1)];

    P_RC = [P_RSW(3,3), P_RSW(3,1);
            P_RSW(1,3), P_RSW(1,1)];

    % ✅ NEW
    P_SW = [P_RSW(2,2), P_RSW(2,3);
            P_RSW(3,2), P_RSW(3,3)];

    %% --- RI ---
    figure(fig_RI);
    h1 = plot(S_err, R_err, 'x', 'Color', color, 'MarkerSize', 8, 'LineWidth', 2);
    handles_RI(end+1) = h1; labels_RI{end+1} = label;
    plot_2D_ellipse(S_err, R_err, P_RI, color);

    %% --- RC ---
    figure(fig_RC);
    h2 = plot(W_err, R_err, 'x', 'Color', color, 'MarkerSize', 8, 'LineWidth', 2);
    handles_RC(end+1) = h2; labels_RC{end+1} = label;
    plot_2D_ellipse(W_err, R_err, P_RC, color);

    %% --- SW (NEW) ---
    figure(fig_SW);
    h3 = plot(S_err, W_err, 'x', 'Color', color, 'MarkerSize', 8, 'LineWidth', 2);
    handles_SW(end+1) = h3; labels_SW{end+1} = label;
    plot_2D_ellipse(S_err, W_err, P_SW, color);

end

%% 5. Finalize
figure(fig_RI);
legend(handles_RI, labels_RI, 'Location','best');
plot(0,0,'k+','HandleVisibility','off');

figure(fig_RC);
legend(handles_RC, labels_RC, 'Location','best');
plot(0,0,'k+','HandleVisibility','off');

figure(fig_SW);
legend(handles_SW, labels_SW, 'Location','best');
plot(0,0,'k+','HandleVisibility','off');

fprintf('\nPlotting Complete.\n');

%% ============================================================
function plot_2D_ellipse(x_center, y_center, P_2x2, color)

    scale = 3;
    [V,D] = eig(P_2x2);

    [d_sorted, idx] = sort(diag(D), 'descend');
    V = V(:,idx);

    a = scale * sqrt(max(0, d_sorted(1)));
    b = scale * sqrt(max(0, d_sorted(2)));

    theta = linspace(0,2*pi,100);
    ellipse_local = [a*cos(theta); b*sin(theta)];
    ellipse_global = V * ellipse_local;

    plot(x_center + ellipse_global(1,:), ...
         y_center + ellipse_global(2,:), ...
         'Color', color, 'LineWidth', 1.5, ...
         'HandleVisibility','off');
end