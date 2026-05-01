function accel_drag_meters_s2 = Atmospheric_Drag_X_Face_Only(r_ECI_meters, v_ECI_meters_s, C_drag)

    mass_kg = SatelliteProperties.SATELLITE_MASS_KG;

    omega_cross_r = [-PhysicsConstants.OMEGA_EARTH_RAD_S*r_ECI_meters(2);...
                      PhysicsConstants.OMEGA_EARTH_RAD_S*r_ECI_meters(1);...
                      0]; 

    vel_relative_meters_s = v_ECI_meters_s - omega_cross_r;

    % v_hat_unit = vel_relative_meters_s/norm(vel_relative_meters_s);
    v_rel_hat_unit = vel_relative_meters_s/norm(vel_relative_meters_s); % Test with not relative velocity
    
    Area_m2 = SatelliteProperties.AREA_X_FACE_M2;

    density_kg_m3 = PhysicsConstants.RHO_0_DRAG_KG_M3 *...
                        exp(-(norm(r_ECI_meters/Units.KILOMETERS) - PhysicsConstants.R_0_DRAG_KM)/PhysicsConstants.H_DRAG_KM); % kg/m^3

    accel_drag_meters_s2 = -(1/2)*C_drag*(Area_m2/mass_kg)*density_kg_m3*norm(vel_relative_meters_s)^2 * v_rel_hat_unit; % (m2/kg)*(kg/m3)*(m/s)*(m/s) = m*m/(m*s^2) = m/s^2
end