function accel_solar_radiation_pressure_meters_s2 = Solar_Radiation_Pressure(r_sat_rel_earth_ECI_meters, r_sun_rel_earth_ECI_meters, sat_is_illuminated) 
    
    % ONLY SOLAR PANEL FORCE!!
    % Vector pointing FROM satellite TO sun
    r_sat_sun_ECI_meters = r_sun_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters;
    dist_sat_sun_meters = norm(r_sat_sun_ECI_meters);
    
    % Unit vectors (pointing towards the sun)
    u_hat = r_sat_sun_ECI_meters / dist_sat_sun_meters; % Incident direction (u_hat in image)
    n_hat = u_hat; % Sun-facing solar panel normal points directly at the sun
    
    phi_inclination_deg = 0; % Angle between surface normal and incident ray
    
    % Heliocentric distance in Astronomical Units (required by the 1/d^2 term)
    AU_METERS = PhysicsConstants.AU_KM * Units.KILOMETERS; 
    d_AU = dist_sat_sun_meters / AU_METERS;

    
    % Calculate reflection coefficients per the image
    nu = (1/3) * SatelliteProperties.SOLAR_CELLS_CD; % Diffuse reflection coefficient, 1/3 * (1 - beta) * gamma (assume gamma = 1)
    mu = (1/2) *  SatelliteProperties.SOLAR_CELLS_CS; % Specular reflection coefficient, (1/2) * Beta * gamma (assume gamma = 1)
    
    % Bidirectional reflectance term
    B_theta = 2 * nu * cosd(phi_inclination_deg) + 4 * mu * cosd(phi_inclination_deg)^2;
    
    % Acceleration calculation (Force / Mass)
    % Added NEGATIVE sign so acceleration pushes AWAY from the sun (opposite to n_hat and s_hat)
    a_srp_meters_s2 = - (PhysicsConstants.SOLAR_PRESSURE_N_M2 / d_AU^2) * ...
                      (SatelliteProperties.AREA_SOLAR_PANEL_M2 / SatelliteProperties.SATELLITE_MASS_KG) * ...
                      ( B_theta * n_hat + (1 - mu) * cosd(phi_inclination_deg)^2 * u_hat );

    % Apply eclipse illumination factor (1 = illuminated, 0 = eclipse)
    accel_solar_radiation_pressure_meters_s2 = a_srp_meters_s2 * sat_is_illuminated;

end