function ref_data = get_reference_data_1Day_Subset()
    persistent ref_data_cached
    if isempty(ref_data_cached)
        fprintf('Loading reference data from .mat files...\n');
        
        % 1. Load and convert to table
        raw_data = load('ref_data/LEO_DATA_Apparent_3Days.mat').LEO_DATA_Apparent;
        Actual_Measurements = array2table(raw_data);
        Actual_Measurements.Properties.VariableNames = {'station_id', 'time_sec_past_epoch', 'apparent_range_km', 'apparent_range_rate_km_s'};
        
        % 2. Filter using table logic (1 day = 86400 seconds)
        day_one_mask = Actual_Measurements.time_sec_past_epoch <= 86400;
        Actual_Measurements = Actual_Measurements(day_one_mask, :);
        
        % 3. Convert the filtered table to a scalar struct for caching
        ref_data_cached.Actual_Measurements = table2struct(Actual_Measurements, 'ToScalar', true);
    end
    ref_data = ref_data_cached;
end