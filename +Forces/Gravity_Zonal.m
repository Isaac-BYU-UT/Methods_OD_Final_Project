function accel_Zonals_km_s2 = Gravity_Zonal(r_km, J2_on, J3_on, J4_on)
    % TODO: Transform to ECEF!
    
    % Set defaults
    if nargin < 2, J2_on = true;  end
    if nargin < 3, J3_on = false; end
    if nargin < 4, J4_on = false; end

    % Constants
    mu = Constants.MU_EARTH_KM3_S2;
    Re = Constants.R_EARTH_KM;
    
    r_mag = sqrt(r_km(1)^2 + r_km(2)^2 + r_km(3)^2);
    phi_sub = r_km(3) / r_mag; % sin(latitude) = z/r -- this needs to be in ECEF!!!!
    
    % Initialize Potential
    U_Zonals = 0;
    
    % J2: Second Degree
    if J2_on
        P2 = 0.5 * (3*phi_sub^2 - 1);
        U_Zonals = U_Zonals - (mu/r_mag) * (Re/r_mag)^2 * Constants.J2_EARTH * P2;
    end
    
    % J3: Third Degree
    if J3_on
        P3 = 0.5 * (5*phi_sub^3 - 3*phi_sub);
        U_Zonals = U_Zonals - (mu/r_mag) * (Re/r_mag)^3 * Constants.J3_EARTH * P3;
    end
    
    % J4: Fourth Degree
    if J4_on
        P4 = (1/8) * (35*phi_sub^4 - 30*phi_sub^2 + 3);
        U_Zonals = U_Zonals - (mu/r_mag) * (Re/r_mag)^4 * Constants.J4_EARTH * P4;
    end

    % Return the gradient vector (Acceleration)
    accel_Zonals_km_s2 = gradient(U_Zonals, r_km); %
end