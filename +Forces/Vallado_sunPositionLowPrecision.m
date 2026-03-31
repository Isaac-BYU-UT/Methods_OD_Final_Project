function [rSun_KM_MOD, rightAscension_rad, declination_rad] = Vallado_sunPositionLowPrecision(jd_utc_days)
% -------------------------------------------------------------------------
% REFERENCE: https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/main/software/matlab/sun.m
% CODE BASED ON VERSION BY VALLADO!
% sunPositionLowPrecision -- adapted from Vallado's code
%
% Computes the geocentric equatorial position of the Sun using Vallado's
% low-precision algorithm (valid ~1950–2050).
%
% INPUT:
%   jd_utc_days          - Julian Date (UTC) [days]
%
% OUTPUTS:
%   rSun_AU_TEME         - Sun position vector [AU] (TEME-like frame)
%   rightAscension_rad   - Right ascension [rad]
%   declination_rad      - Declination [rad]
%
% NOTES:
%   - Accuracy ~0.01 deg
%   - Internally uses degrees where Vallado specifies, then converts to rad
% -------------------------------------------------------------------------

    %% -------------------- Constants --------------------
    TWO_PI_rad   = 2.0 * pi;
    DEG2RAD      = pi / 180.0;

    %% -------------------- Time Variables --------------------
    % Julian centuries since J2000 (UT1 approximation)
    T_UT1_centuries = (jd_utc_days - 2451545.0) / 36525.0;

    % For this model, assume TDB ≈ UT1 !!!!!!!!!
    T_TDB_centuries = T_UT1_centuries;

    %% -------------------- Solar Angles (Degrees) --------------------
    % Mean longitude of the Sun [deg]
    meanLongitude_deg = 280.460 + 36000.771285 * T_UT1_centuries;
    meanLongitude_deg = mod(meanLongitude_deg, 360.0);

    % Mean anomaly of the Sun [rad]
    meanAnomaly_deg = 357.528 + 35999.0509575 * T_TDB_centuries;
    meanAnomaly_rad = mod(meanAnomaly_deg * DEG2RAD, TWO_PI_rad);

    if meanAnomaly_rad < 0.0
        meanAnomaly_rad = meanAnomaly_rad + TWO_PI_rad;
    end

    % Ecliptic longitude [deg]
    eclipticLongitude_deg = meanLongitude_deg ...
        + 1.914666471 * sin(meanAnomaly_rad) ...
        + 0.019994643 * sin(2.0 * meanAnomaly_rad);
    eclipticLongitude_deg = mod(eclipticLongitude_deg, 360.0);

    % Mean obliquity of the ecliptic [deg]
    obliquity_deg = 23.439291 - 0.0130042 * T_TDB_centuries;

    %% -------------------- Convert to Radians --------------------
    eclipticLongitude_rad = eclipticLongitude_deg * DEG2RAD;
    obliquity_rad         = obliquity_deg * DEG2RAD;

    %% -------------------- Sun-Earth (Magnitude) --------------------
    % Distance from Earth to Sun [AU]
    rSun_AU = 1.000140612 ...
        - 0.016708617 * cos(meanAnomaly_rad) ...
        - 0.000139589 * cos(2.0 * meanAnomaly_rad);

    %% -------------------- Position Vector (TEME-like) --------------------
    rSun_AU_MOD = zeros(3,1);

    rSun_AU_MOD(1) = rSun_AU * cos(eclipticLongitude_rad);
    rSun_AU_MOD(2) = rSun_AU * cos(obliquity_rad) * sin(eclipticLongitude_rad);
    rSun_AU_MOD(3) = rSun_AU * sin(obliquity_rad) * sin(eclipticLongitude_rad);

    rSun_KM_MOD = rSun_AU_MOD * Constants.AU_KM;

    %% -------------------- Right Ascension --------------------
    rightAscension_rad = atan( cos(obliquity_rad) * tan(eclipticLongitude_rad) );

    % Ensure correct quadrant
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