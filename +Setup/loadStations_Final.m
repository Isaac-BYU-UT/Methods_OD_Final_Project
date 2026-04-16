function Stations = loadStations_Final()

    raw_station_positions_ECEF_m = [
        -6143584,  1364250,  1033743;
        1907295,  6030810,  -817119;
        2390310, -5564341,  1994578];

    % Station 1
    Stations(1).name = 'Atoll';
    Stations(1).position_ECEF_meters = raw_station_positions_ECEF_m(1, :)';
    Stations(1).sigma_range_meters = 10; % Standard deviation of range measurement in meters
    Stations(1).sigma_range_rate_meters_per_sec = 0.5 * Units.MILIMETERS; % Standard deviation of range rate measurement in meters per second
    Stations(1).Covariance = diag([Stations(1).sigma_range_meters^2, Stations(1).sigma_range_rate_meters_per_sec^2]);

    % Station 2
    Stations(2).name = 'Diego Garcia';
    Stations(2).position_ECEF_meters = raw_station_positions_ECEF_m(2, :)';
    Stations(2).sigma_range_meters = 5; % Standard deviation of range measurement in meters
    Stations(2).sigma_range_rate_meters_per_sec = 1.0 * Units.MILIMETERS; % Standard deviation of range rate measurement in meters per second
    Stations(2).Covariance = diag([Stations(2).sigma_range_meters^2, Stations(2).sigma_range_rate_meters_per_sec^2]);

    % Station 3
    Stations(3).name = 'Arecibo';
    Stations(3).position_ECEF_meters = raw_station_positions_ECEF_m(3, :)';
    Stations(3).sigma_range_meters = 10 ; % Standard deviation of range measurement in meters
    Stations(3).sigma_range_rate_meters_per_sec = 0.5 * Units.MILIMETERS; % Standard deviation of range rate measurement in meters per second
    Stations(3).Covariance = diag([Stations(3).sigma_range_meters^2, Stations(3).sigma_range_rate_meters_per_sec^2]);

end