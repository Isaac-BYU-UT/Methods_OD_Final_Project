function ref_data = get_reference_data_All6Days()

    persistent ref_data_cached

    if isempty(ref_data_cached)
        fprintf('Loading reference data from .mat files...\n');

        Actual_Measurements_A = load('ref_data/LEO_DATA_Apparent_3Days.mat').LEO_DATA_Apparent;

        Actual_Measurements_A = array2table(Actual_Measurements_A); % Convert to table for easier handling
        Actual_Measurements_A.Properties.VariableNames = {'station_id', 'time_sec_past_epoch', 'apparent_range_km', 'apparent_range_rate_km_s'}; % Make the variable names more intuitive

        Actual_Measurements_B = load('ref_data/LEO_DATA_Apparent_Days4-6.mat').LEO_DATA_Apparent;
        
        Actual_Measurements_B = array2table(Actual_Measurements_B); % Convert to table for easier handling
        Actual_Measurements_B.Properties.VariableNames = {'station_id', 'time_sec_past_epoch', 'apparent_range_km', 'apparent_range_rate_km_s'}; % Make the variable names more intuitive

        Actual_Measurements_A = struct('Actual_Measurements', table2struct(Actual_Measurements_A,'ToScalar',true));
        Actual_Measurements_B = struct('Actual_Measurements', table2struct(Actual_Measurements_B,'ToScalar',true));

        Actual_Measurements_Combined.station_id = [Actual_Measurements_A.Actual_Measurements.station_id; Actual_Measurements_B.Actual_Measurements.station_id];
        Actual_Measurements_Combined.time_sec_past_epoch = [Actual_Measurements_A.Actual_Measurements.time_sec_past_epoch; ...
                                                            Actual_Measurements_B.Actual_Measurements.time_sec_past_epoch + 3 * Units.SEC_IN_SOLAR_DAY];
        Actual_Measurements_Combined.apparent_range_km = [Actual_Measurements_A.Actual_Measurements.apparent_range_km; Actual_Measurements_B.Actual_Measurements.apparent_range_km];
        Actual_Measurements_Combined.apparent_range_rate_km_s = [Actual_Measurements_A.Actual_Measurements.apparent_range_rate_km_s; Actual_Measurements_B.Actual_Measurements.apparent_range_rate_km_s];

        Final_Meas.Actual_Measurements = Actual_Measurements_Combined;

        ref_data_cached = Final_Meas;

    end

    ref_data = ref_data_cached;
    
end