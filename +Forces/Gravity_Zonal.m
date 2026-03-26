function accel_Zonals_km_s2 = Gravity_Zonal(r_km, J2_inc, J3_inc, J4_inc)
    mu = Constants.MU_EARTH_KM3_S2;
    R_Earth = Constants.R_EARTH_KM;
    J2 = Constants.J2_EARTH;
    J3 = Constants.J3_EARTH;
    J4 = Constants.J4_EARTH;

    if nargin < 2, J2_inc = true; end % Default is juts J2
    if nargin < 3, J3_inc = false; end
    if nargin < 4, J4_inc = false; end

    U_Zonals = 0;
    if J2_inc
        U_Zonals = U_Zonals + Potential_J2(r_km, mu, R_Earth, J2);
    end
    if J3_inc
        U_Zonals = U_Zonals + Potential_J3(r_km, mu, R_Earth, J3);
    end
    if J4_inc
        U_Zonals = U_Zonals + Potential_J4(r_km, mu, R_Earth, J4);
    end

    accel_Zonals_km_s2 = gradient(U_Zonals, r_km);
end


function U = Potential_J2(r, mu, Re, J2) % Vallado pg. 594
    r_mag = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    phi_sub = r(3) / r_mag; % sin of geocentric latitude
    U = (mu / r_mag) * (Re/r_mag)^2 * J2 * (1/2) * (3*phi_sub^2 - 1);
end

function U = Potential_J3(r, mu, Re, J3) % Vallado pg. 594
    r_mag = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    phi_sub = r(3) / r_mag;
    U = (mu / r_mag) * (Re/r_mag)^3 * J3 * (1/2) * (5*phi_sub^3 - 3*phi_sub);
end

function U = Potential_J4(r, mu, Re, J4) % Vallado pg. 594
    
    r_mag = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    phi_sub = r(3) / r_mag; % sin of geocentric latitude (z/r)
    
    U = -(mu / r_mag) * (Re / r_mag)^4 * J4 * (1/8) * (35*phi_sub^4 - 30*phi_sub^2 + 3);
end