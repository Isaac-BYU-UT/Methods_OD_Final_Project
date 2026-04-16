function Measurement_Table = make_measurement_table(t_obs, obs_measurements, Y_computed)
    % Add computed values and pre-fit residuals to the Measurement_Table

    Measurement_Table.time_sec_past_epoch = t_obs;
    Measurement_Table.apparent_range_meters = obs_measurements.apparent_range_km * Units.KILOMETERS;
    Measurement_Table.computed_range_meters = Y_computed(:,1);
    Measurement_Table.apparent_range_rate_meters_s = obs_measurements.apparent_range_rate_km_s * Units.KILOMETERS;
    Measurement_Table.computed_range_rate_meters_s = Y_computed(:,2);
    Measurement_Table.station_id = obs_measurements.station_id;
end