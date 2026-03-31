function isIlluminated = Vallado_sunLOS(r_sat_km, r_sun_km)
% -------------------------------------------------------------------------
% Determines if a satellite has line-of-sight to the Sun (i.e., is illuminated)
% using Vallado Algorithm 35 (SIGHT) -- Cylindrical Model of Shadow
%
% Inputs:
%   r_sat_km : [3x1] satellite position vector in ECI (km)
%   r_sun_km : [3x1] Sun position vector in ECI (km)
%
% Output:
%   isIlluminated : boolean (true = sunlit, false = in Earth's shadow)
%
% -------------------------------------------------------------------------

    % Earth's radius (km)
    R_earth = Constants.R_EARTH_KM;

    % Precompute dot products
    r1 = r_sat_km(:);
    r2 = r_sun_km(:);

    r1_dot_r1 = dot(r1, r1);
    r2_dot_r2 = dot(r2, r2);
    r1_dot_r2 = dot(r1, r2);

    % Compute tau_min
    tau_min = (r1_dot_r1 - r1_dot_r2) / ...
              (r1_dot_r1 + r2_dot_r2 - 2 * r1_dot_r2);

    % disp(tau_min);

    % Default: assume NO line of sight (i.e., eclipse)
    isIlluminated = false;

    % Case 1: closest approach not between the two points
    if (tau_min < 0) || (tau_min > 1)
        isIlluminated = true;
        return;
    end

    % Case 2: compute closest approach distance squared
    c_squared = (1 - tau_min) * r1_dot_r1 + tau_min * r1_dot_r2;

    % disp(c2/R_earth^2);

    % If closest distance is outside Earth → LOS exists
    if c_squared >= R_earth^2
        isIlluminated = true;
    end

end