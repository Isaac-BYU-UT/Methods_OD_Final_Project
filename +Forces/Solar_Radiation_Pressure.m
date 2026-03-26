function accel_solar_radiation_pressure_km_s2 = Solar_Radiation_Pressure(r_sat_rel_earth_ECI_km, r_sun_rel_earth_ECI_km)
    
    r_sat_sun = r_sun_rel_earth_ECI_km - r_sat_rel_earth_ECI_km;
    a_rad_sun = -1* (Constants.SOLAR_PRESSURE_N_M2 * Constants.SRP_AREA_M2 * Constants.C_Reflectivity)...
                /Constants.SATTELITE_MASS_KG  * ...
                 (r_sat_sun/norm(r_sat_sun)^3);

    %TODO: Compute projected area
    % TODO: Compute reflectivity coefficient
    % TODO: Check Earth's Shadow (is satellite in shadow or not?)

    accel_solar_radiation_pressure_km_s2 = a_rad_sun;
end