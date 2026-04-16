function ref_data = get_reference_data_3Days()

    persistent ref_data_cached

    if isempty(ref_data_cached)
        fprintf('Loading reference data from .mat files...\n');

        Actual_Measurements = load('ref_data/LEO_DATA_Apparent_3Days.mat').LEO_DATA_Apparent;
        Actual_Measurements = array2table(Actual_Measurements); % Convert to table for easier handling

        % Convert units from km to m and km/s to m/s
        Actual_Measurements(:, 3) = Actual_Measurements(:, 3) * Units.KILOMETERS; % Scale range to meters
        Actual_Measurements(:, 4) = Actual_Measurements(:, 4) * Units.KILOMETERS; % Scale range rate to meters/second
        Actual_Measurements.Properties.VariableNames = {'station_id', 'time_sec_past_epoch', 'apparent_range_meters', 'apparent_range_rate_meters_per_second'}; % Make the variable names more intuitive

        ref_data_cached = struct('Actual_Measurements', table2struct(Actual_Measurements,'ToScalar',true));

    end

    ref_data = ref_data_cached;
    
end