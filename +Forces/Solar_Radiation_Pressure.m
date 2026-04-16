function accel_solar_radiation_pressure_meters_s2 = Solar_Radiation_Pressure(r_sat_rel_earth_ECI_meters, r_sun_rel_earth_ECI_meters,sat_is_illuminated)
    
% ONLY SOLAR PANEL FORCE!!
r_sat_sun_ECI_meters = r_sun_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters;


    phi_inclination_deg = 0; % Angle between suface normal and incoming radiation, 0 deg for solar panel
    n_hat = r_sat_sun_ECI_meters/norm(r_sat_sun_ECI_meters);
    s_hat = r_sat_sun_ECI_meters/norm(r_sat_sun_ECI_meters);
    
    a_srp_meters_s2 = (Constants.SOLAR_PRESSURE_N_M2 * Constants.AREA_SOLAR_PANEL_M2 * cosd(phi_inclination_deg) / Constants.SATTELITE_MASS_KG) *...
              (...
                2*(Constants.SOLAR_CELLS_CD / 3 + Constants.SOLAR_CELLS_CS*cosd(phi_inclination_deg)) * n_hat + ...
                (1- Constants.SOLAR_CELLS_CS)*s_hat ...
                );

    accel_solar_radiation_pressure_meters_s2 = a_srp_meters_s2 * sat_is_illuminated; % Will switch to 0 if not illuminated

end