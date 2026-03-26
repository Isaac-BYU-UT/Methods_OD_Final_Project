function ref_data = get_reference_data()

    persistent ref_data_cached

    if isempty(ref_data_cached)
        fprintf('Loading reference data from .mat files...\n');
        % Load the reference data from the .mat file
        A_t0_ref = load('ref_data/A_t0.mat').A;
        H_tilde_t0_ref = load('ref_data/H_Tilde_t0.mat').H_TILDA;
        STM_21600_0_ref = load('ref_data/Phi_21600_0.mat').PHI_t_120;
        Light_Time_Truth_ref = load('ref_data/lighttime_truth.mat').lighttime_truth;
        Light_Time_Truth_ref = array2table(Light_Time_Truth_ref);
        Light_Time_Truth_ref.Properties.VariableNames = {'r_x_km', 'r_y_km', 'r_z_km', 'v_x_km_s', 'v_y_km_s', 'v_z_km_s', 'time_since_epoch_sec'}; % Make the variable names more intuitive


        Actual_Measurements = load('ref_data/LEO_DATA_Apparent.mat').LEO_DATA_Apparent;
        Actual_Measurements = array2table(Actual_Measurements); % Convert to table for easier handling
        Actual_Measurements.Properties.VariableNames = {'station_id', 'time_sec_past_epoch', 'apparent_range_km', 'apparent_range_rate_km_s'}; % Make the variable names more intuitive



        ref_data_cached = struct('A_t0_ref', A_t0_ref, 'H_tilde_t0_ref', H_tilde_t0_ref, 'STM_21600_0_ref', STM_21600_0_ref, 'Light_Time_Truth_ref', Light_Time_Truth_ref, 'Actual_Measurements', Actual_Measurements);

    end

    ref_data = ref_data_cached;
    
end