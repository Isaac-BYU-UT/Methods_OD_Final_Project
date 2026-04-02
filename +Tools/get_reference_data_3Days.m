function ref_data = get_reference_data_3Days()

    persistent ref_data_cached

    if isempty(ref_data_cached)
        fprintf('Loading reference data from .mat files...\n');

        Actual_Measurements = load('ref_data/LEO_DATA_Apparent_3Days.mat').LEO_DATA_Apparent;
        Actual_Measurements = array2table(Actual_Measurements); % Convert to table for easier handling
        Actual_Measurements.Properties.VariableNames = {'station_id', 'time_sec_past_epoch', 'apparent_range_km', 'apparent_range_rate_km_s'}; % Make the variable names more intuitive

        ref_data_cached = struct('Actual_Measurements', Actual_Measurements);

    end

    ref_data = ref_data_cached;
    
end