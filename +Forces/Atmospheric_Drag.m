function accel_drag_km_s2 = Atmospheric_Drag(r_km, v_km_s, C_drag)

    mass_kg = Constants.SATTELITE_MASS_KG; % kg
    Area_m2 = Constants.SATELLITE_AREA_M2; % m^2

    density_kg_m3 = Constants.RHO_0_DRAG_KG_M3 * exp(-(norm(r_km) - Constants.R_0_DRAG_KM)/Constants.H_DRAG_KM); % kg/m^3
    vel_relative_km_s = v_km_s + [Constants.OMEGA_EARTH_RAD_S*r_km(2); -Constants.OMEGA_EARTH_RAD_S*r_km(1); 0]; % Account for earth's rotational velocity, TODO; Could make this more accurate with LOD?
    accel_drag_km_s2 = -(1/2)*C_drag*(Area_m2/mass_kg)*density_kg_m3*vel_relative_km_s*norm(vel_relative_km_s) * 1000; % (m2/kg)*(kg/m3)*(km/s)*(km/s) = km*km/(m*s^2) * (1000m/km) = km/s^2

    % TODO: Make Area projected based on velocity.
    

end
