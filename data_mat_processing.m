clear; clc;

results_dir = 'Results/4_23_2026/';

% Initialize all cases as NaN (default)
sorensen_pos_caseA = NaN(3,1); sorensen_poscov_caseA = NaN(6,6);
sorensen_pos_caseB = NaN(3,1); sorensen_poscov_caseB = NaN(6,6);
sorensen_pos_caseC = NaN(3,1); sorensen_poscov_caseC = NaN(6,6);
sorensen_pos_caseD = NaN(3,1); sorensen_poscov_caseD = NaN(6,6);
sorensen_pos_caseE = NaN(3,1); sorensen_poscov_caseE = NaN(6,6);
sorensen_pos_caseF = NaN(3,1); sorensen_poscov_caseF = NaN(6,6);
sorensen_pos_caseG = NaN(3,1); sorensen_poscov_caseG = NaN(6,6);

files = dir(fullfile(results_dir, '*.mat'));

for i = 1:length(files)
    filename = files(i).name;
    filepath = fullfile(results_dir, filename);
    
    data = load(filepath); % load everything safely
    
    % Default values if fields are missing
    if isfield(data, 'r_final_ECI_meters')
        r = data.r_final_ECI_meters;
    else
        r = NaN(3,1);
    end
    
    if isfield(data, ' P_cov_final')
        P = data. P_cov_final;
    else
        P = NaN(6,6);
    end

    % Case (a): Range only, all stations
    if contains(filename, '_A_')
        sorensen_pos_caseA = r;
        sorensen_poscov_caseA = P;
    end
    
    % Case (b): Range + Range-rate, all stations
    if contains(filename, '_B_')
        sorensen_pos_caseB = r;
        sorensen_poscov_caseB = P;
    end
    
    % Case (c): Kwajalein (Atoll)
    if contains(filename, '_C_')
        sorensen_pos_caseC = r;
        sorensen_poscov_caseC = P;
    end
    
    % Case (d): Diego Garcia
    if contains(filename, '_D_')
        sorensen_pos_caseD = r;
        sorensen_poscov_caseD = P;
    end
    
    % Case (e): Arecibo
    if contains(filename, '_E_')
        sorensen_pos_caseE = r;
        sorensen_poscov_caseE = P;
    end

    % Case (f): Long Arc
    if contains(filename, '_F_')
        sorensen_pos_caseG = r;
        sorensen_poscov_caseG = P;
    end

    % Case (g): Short Arc (All station, all data)
    if contains(filename, '_G_')
        sorensen_pos_caseG = r;
        sorensen_poscov_caseG = P;
    end
end