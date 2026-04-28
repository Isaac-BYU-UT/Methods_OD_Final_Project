UTC.year = 2018; UTC.month = 2; UTC.day = 1; UTC.hour = 5; UTC.minute = 0; UTC.seconds = 0.0; % Gregorian Date (UTC) 1 Feb 2018, 05:00:00 UTC.
UTC_date_time_epoch_t0 = datetime(UTC.year, UTC.month, UTC.day,UTC.hour, UTC.minute, UTC.seconds, 'TimeZone','UTC');

time_struct = Tools.ComputeTimeSystems(UTC_date_time_epoch_t0);

% Example: Position of the Sun and Moon in ECI (km)
jd = time_struct.jd_UTC_days;
[sunPos, ~] = planetEphemeris(jd, 'Earth', 'Sun');
[moonPos, ~] = planetEphemeris(jd, 'Earth', 'Moon');

% Result is in kilometers
fprintf('Sun ECI Pos: %f %f %f km\n', sunPos);
fprintf('Moon ECI Pos: %f %f %f km\n', moonPos);
