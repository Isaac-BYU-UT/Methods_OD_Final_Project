function accel_luni_solar_meters_s2 = Luni_Solar_Pertubations(r_sat_rel_earth_ECI_meters, r_sun_rel_earth_ECI_meters, r_moon_rel_earth_ECI_meters)
    
    r_sat_moon_ECI_meters = r_moon_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters;
    r_sat_sun_ECI_meters = r_sun_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters;

    a_grav_moon_m_s2 = (PhysicsConstants.MU_MOON_KM3_S2 * Units.KILOMETERS^3)  * ...
                        (r_sat_moon_ECI_meters/norm(r_sat_moon_ECI_meters)^3 - r_moon_rel_earth_ECI_meters/norm(r_moon_rel_earth_ECI_meters)^3);
    
    a_grav_sun_m_s2 = (PhysicsConstants.MU_SUN_KM3_S2 * Units.KILOMETERS^3)  * ...
                    (r_sat_sun_ECI_meters/norm(r_sat_sun_ECI_meters)^3 - r_sun_rel_earth_ECI_meters/norm(r_sun_rel_earth_ECI_meters)^3);

    accel_luni_solar_meters_s2 = a_grav_moon_m_s2 + a_grav_sun_m_s2;
end