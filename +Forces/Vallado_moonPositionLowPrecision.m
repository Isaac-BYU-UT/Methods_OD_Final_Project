function [rMoon_km_MOD, rightAscension_rad, declination_rad] = Vallado_moonPositionLowPrecision(jd_utc_days)
% -------------------------------------------------------------------------
% REFERENCE: https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/main/software/matlab/moon.m 
% CODE BASED ON VERSION BY VALLADO!
% moonPositionLowPrecision
%
% Computes the geocentric equatorial position of the Moon using Vallado's
% low-precision analytical model.
%
% INPUT:
%   jd_utc_days          - Julian Date [days]
%
% OUTPUTS:
%   rMoon_km_ECI         - Moon position vector [km] (ECI frame)
%   rightAscension_rad   - Right ascension [rad]
%   declination_rad      - Declination [rad]
%
% NOTES:
%   - Based on Vallado (Alg 31)
%   - Low-precision (~0.3 deg typical accuracy)
%   - Internally uses degree-based series expansions
% -------------------------------------------------------------------------

    %% -------------------- Constants --------------------
    TWO_PI_rad = 2.0 * pi;
    DEG2RAD    = pi / 180.0;

    % Load astronomical constants (expects Earth radius, etc.)
    re = Constants.R_EARTH_KM;

    %% -------------------- Time --------------------
    % Julian centuries since J2000 (TDB approximation!!!)
    T_TDB_centuries = (jd_utc_days - 2451545.0) / 36525.0;

    %% -------------------- Ecliptic Longitude [deg] --------------------
    eclipticLongitude_deg = 218.32 + 481267.8813 * T_TDB_centuries ...
        + 6.29 * sin((134.9 + 477198.85 * T_TDB_centuries) * DEG2RAD) ...
        - 1.27 * sin((259.2 - 413335.38 * T_TDB_centuries) * DEG2RAD) ...
        + 0.66 * sin((235.7 + 890534.23 * T_TDB_centuries) * DEG2RAD) ...
        + 0.21 * sin((269.9 + 954397.70 * T_TDB_centuries) * DEG2RAD) ...
        - 0.19 * sin((357.5 + 35999.05  * T_TDB_centuries) * DEG2RAD) ...
        - 0.11 * sin((186.6 + 966404.05 * T_TDB_centuries) * DEG2RAD);

    %% -------------------- Ecliptic Latitude [deg] --------------------
    eclipticLatitude_deg = ...
          5.13 * sin(( 93.3 + 483202.03 * T_TDB_centuries) * DEG2RAD) ...
        + 0.28 * sin((228.2 + 960400.87 * T_TDB_centuries) * DEG2RAD) ...
        - 0.28 * sin((318.3 +   6003.18 * T_TDB_centuries) * DEG2RAD) ...
        - 0.17 * sin((217.6 - 407332.20 * T_TDB_centuries) * DEG2RAD);

    %% -------------------- Horizontal Parallax [deg] --------------------
    horizontalParallax_deg = 0.9508 ...
        + 0.0518 * cos((134.9 + 477198.85 * T_TDB_centuries) * DEG2RAD) ...
        + 0.0095 * cos((259.2 - 413335.38 * T_TDB_centuries) * DEG2RAD) ...
        + 0.0078 * cos((235.7 + 890534.23 * T_TDB_centuries) * DEG2RAD) ...
        + 0.0028 * cos((269.9 + 954397.70 * T_TDB_centuries) * DEG2RAD);

    %% -------------------- Convert to Radians --------------------
    eclipticLongitude_rad = mod(eclipticLongitude_deg * DEG2RAD, TWO_PI_rad);
    eclipticLatitude_rad  = mod(eclipticLatitude_deg  * DEG2RAD, TWO_PI_rad);
    horizontalParallax_rad = mod(horizontalParallax_deg * DEG2RAD, TWO_PI_rad);

    %% -------------------- Obliquity of the Ecliptic --------------------
    obliquity_deg = 23.439291 - 0.0130042 * T_TDB_centuries;
    obliquity_rad = obliquity_deg * DEG2RAD;

    %% -------------------- Direction Cosines --------------------
    % Unit vector from Earth to Moon in ECI frame
    l_hat = cos(eclipticLatitude_rad) * cos(eclipticLongitude_rad);

    m_hat = cos(obliquity_rad) * cos(eclipticLatitude_rad) * sin(eclipticLongitude_rad) ...
          - sin(obliquity_rad) * sin(eclipticLatitude_rad);

    n_hat = sin(obliquity_rad) * cos(eclipticLatitude_rad) * sin(eclipticLongitude_rad) ...
          + cos(obliquity_rad) * sin(eclipticLatitude_rad);

    %% -------------------- Distance (km) --------------------
    % Vallado formulation: distance = Earth radius / sin(horizontal parallax)
    rMoon_km = re / sin(horizontalParallax_rad);

    %% -------------------- Position Vector --------------------
    rMoon_km_MOD = zeros(3,1);

    rMoon_km_MOD(1) = rMoon_km * l_hat;
    rMoon_km_MOD(2) = rMoon_km * m_hat;
    rMoon_km_MOD(3) = rMoon_km * n_hat;

    %% -------------------- Right Ascension & Declination --------------------
    rightAscension_rad = atan2(m_hat, l_hat);
    declination_rad    = asin(n_hat);

end