function [rSun_KM_MOD, rightAscension_rad, declination_rad] = Vallado_sunAlmanacLowPrecision(jd_utc_days)
% -------------------------------------------------------------------------
% Reference: https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/main/software/matlab/sunalmanac.m
% CODE BASED ON VERSION BY VALLADO!
% sunAlmanacLowPrecision
%
% Computes the geocentric equatorial position of the Sun using the
% Almanac-based low-precision model (Vallado, ~2002 version).
%
% INPUT:
%   jd_utc_days          - Julian Date [days]
%
% OUTPUTS:
%   rSun_AU_ECI          - Sun position vector [AU] (ECI-like frame)
%   rightAscension_rad   - Right ascension [rad]
%   declination_rad      - Declination [rad]
%
% NOTES:
%   - Valid approximately from 1950 to 2050
%   - Accuracy ~0.01 degrees
%   - Uses degree-based almanac formulation internally
% -------------------------------------------------------------------------

    %% -------------------- Constants --------------------
    TWO_PI_rad = 2.0 * pi;
    DEG2RAD    = pi / 180.0;

    %% -------------------- Time --------------------
    % Julian centuries since J2000 (UT1 approximation)
    T_UT1_centuries = (jd_utc_days - 2451545.0) / 36525.0;

    % Approximate TDB ≈ UT1 for this low-precision model
    T_TDB_centuries = T_UT1_centuries;

    %% -------------------- Mean Solar Quantities (Degrees) --------------------
    % Mean longitude of the Sun [deg]
    meanLongitude_deg = 280.460 + 36000.771285 * T_UT1_centuries;
    meanLongitude_deg = mod(meanLongitude_deg, 360.0);

    % Mean anomaly [rad]
    meanAnomaly_deg = 357.528 + 35999.050957 * T_TDB_centuries;
    meanAnomaly_rad = mod(meanAnomaly_deg * DEG2RAD, TWO_PI_rad);

    if meanAnomaly_rad < 0.0
        meanAnomaly_rad = meanAnomaly_rad + TWO_PI_rad;
    end

    % Ecliptic longitude [deg]
    eclipticLongitude_deg = meanLongitude_deg ...
        + 1.915 * sin(meanAnomaly_rad) ...
        + 0.020 * sin(2.0 * meanAnomaly_rad);
    eclipticLongitude_deg = mod(eclipticLongitude_deg, 360.0);

    % Mean obliquity of the ecliptic [deg]
    obliquity_deg = 23.439 - 0.01461 * T_TDB_centuries;

    %% -------------------- Convert to Radians --------------------
    eclipticLongitude_rad = eclipticLongitude_deg * DEG2RAD;
    obliquity_rad         = obliquity_deg * DEG2RAD;

    %% -------------------- Sun-Earth Distance --------------------
    % Distance from Earth to Sun [AU]
    rSun_AU = 1.00014 ...
        - 0.01671 * cos(meanAnomaly_rad) ...
        - 0.00014 * cos(2.0 * meanAnomaly_rad);

    %% -------------------- Position Vector --------------------
    rSun_AU_MOD = zeros(3,1);

    rSun_AU_MOD(1) = rSun_AU * cos(eclipticLongitude_rad);
    rSun_AU_MOD(2) = rSun_AU * cos(obliquity_rad) * sin(eclipticLongitude_rad);
    rSun_AU_MOD(3) = rSun_AU * sin(obliquity_rad) * sin(eclipticLongitude_rad);

    rSun_KM_MOD = rSun_AU_MOD * Constants.AU_KM;

    %% -------------------- Right Ascension --------------------
    rightAscension_rad = atan( cos(obliquity_rad) * tan(eclipticLongitude_rad) );

    % Quadrant correction
    if eclipticLongitude_rad < 0.0
        eclipticLongitude_rad = eclipticLongitude_rad + TWO_PI_rad;
    end

    if abs(eclipticLongitude_rad - rightAscension_rad) > (pi / 2)
        rightAscension_rad = rightAscension_rad + ...
            (pi/2) * round((eclipticLongitude_rad - rightAscension_rad) / (pi/2));
    end

    %% -------------------- Declination --------------------
    declination_rad = asin( sin(obliquity_rad) * sin(eclipticLongitude_rad) );

end