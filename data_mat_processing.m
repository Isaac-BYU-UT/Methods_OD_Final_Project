clear; clc;

results_dir = 'Results/4_2_2026/';

% Initialize all cases as NaN (default)
sorensen_pos_caseA = NaN(3,1); sorensen_poscov_caseA = NaN(7,7);
sorensen_pos_caseB = NaN(3,1); sorensen_poscov_caseB = NaN(7,7);
sorensen_pos_caseC = NaN(3,1); sorensen_poscov_caseC = NaN(7,7);
sorensen_pos_caseD = NaN(3,1); sorensen_poscov_caseD = NaN(7,7);
sorensen_pos_caseE = NaN(3,1); sorensen_poscov_caseE = NaN(7,7);
sorensen_pos_caseF = NaN(3,1); sorensen_poscov_caseF = NaN(7,7);
sorensen_pos_caseG = NaN(3,1); sorensen_poscov_caseG = NaN(7,7);

files = dir(fullfile(results_dir, '*.mat'));

for i = 1:length(files)
    filename = files(i).name;
    filepath = fullfile(results_dir, filename);
    
    data = load(filepath); % load everything safely
    
    % Default values if fields are missing
    if isfield(data, 'r_sat_t_dV1_ECI_km')
        r = data.r_sat_t_dV1_ECI_km;
    else
        r = NaN(3,1);
    end
    
    if isfield(data, 'P_bar_cov_dV1')
        P = data.P_bar_cov_dV1;
    else
        P = NaN(7,7);
    end

    % Case (a): Range only, all stations
    if contains(filename, 'All_station_Range_data')
        sorensen_pos_caseA = r;
        sorensen_poscov_caseA = P;
    end
    
    % Case (b): Range + Range-rate, all stations
    if contains(filename, 'All_station_Range_Rate_data')
        sorensen_pos_caseB = r;
        sorensen_poscov_caseB = P;
    end
    
    % Case (c): Kwajalein (Atoll)
    if contains(filename, 'Atoll')
        sorensen_pos_caseC = r;
        sorensen_poscov_caseC = P;
    end
    
    % Case (d): Diego Garcia
    if contains(filename, 'Diego_Garcia')
        sorensen_pos_caseD = r;
        sorensen_poscov_caseD = P;
    end
    
    % Case (e): Arecibo
    if contains(filename, 'Arecibo')
        sorensen_pos_caseE = r;
        sorensen_poscov_caseE = P;
    end

    % Case (g): Short Arc (All station, all data)
    if contains(filename, 'All_station_All_data')
        sorensen_pos_caseG = r;
        sorensen_poscov_caseG = P;
    end
end

% Case (f): Long Arc (not implemented yet → keep NaNs)
sorensen_pos_caseF = NaN(3,1);
sorensen_poscov_caseF = NaN(7,7);