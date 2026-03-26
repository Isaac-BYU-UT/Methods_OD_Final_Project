function accel_luni_solar_km_s2 = Luni_Solar_Pertubations(r_sat_rel_earth_ECI_km, r_sun_rel_earth_ECI_km, r_moon_rel_earth_ECI_km)
    
    r_sat_moon = r_moon_rel_earth_ECI_km - r_sat_rel_earth_ECI_km;
    r_sat_sun = r_sun_rel_earth_ECI_km - r_sat_rel_earth_ECI_km;
    a_grav_moon = Constants.MU_MOON_KM3_S2  * (r_sat_moon/norm(r_sat_moon)^3 - r_moon_rel_earth_ECI_km/norm(r_moon_rel_earth_ECI_km)^3);
    a_grav_sun = Constants.MU_SUN_KM3_S2  * (r_sat_sun/norm(r_sat_sun)^3 - r_sun_rel_earth_ECI_km/norm(r_sun_rel_earth_ECI_km)^3);

    accel_luni_solar_km_s2 = a_grav_moon + a_grav_sun;
end