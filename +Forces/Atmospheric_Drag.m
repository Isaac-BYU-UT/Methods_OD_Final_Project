function accel_drag_meters_s2 = Atmospheric_Drag(r_ECI_meters, v_ECI_meters_s, C_drag, r_sun_rel_earth_ECI_meters)

    % TODO: Unit test this!

    mass_kg = SatelliteProperties.SATELLITE_MASS_KG;

    % Account for earth's rotational velocity, TODO; Could make this more accurate with LOD
    omega_cross_r = [PhysicsConstants.OMEGA_EARTH_RAD_S*r_ECI_meters(2);...
                    -PhysicsConstants.OMEGA_EARTH_RAD_S*r_ECI_meters(1);...
                     0]; 

    vel_relative_meters_s = v_ECI_meters_s - omega_cross_r;

    % v_hat_unit = vel_relative_meters_s/norm(vel_relative_meters_s);
    v_hat_unit = vel_relative_meters_s/norm(vel_relative_meters_s); % Test with not relative velocity
    
    angular_momentum_vec_ECI = Tools.crossProductEquivalent(r_ECI_meters) * v_ECI_meters_s;
    r_sat_to_sun_vec_ECI = r_sun_rel_earth_ECI_meters - r_ECI_meters;

    z_hat_unit_vec_ECI = r_ECI_meters/norm(r_ECI_meters);
    y_hat_unit_vec_ECI = -1*angular_momentum_vec_ECI/norm(angular_momentum_vec_ECI);
    x_hat_unit_vec_ECI = Tools.crossProductEquivalent(y_hat_unit_vec_ECI) * z_hat_unit_vec_ECI;
    solar_panel_unit_vec = r_sat_to_sun_vec_ECI/norm(r_sat_to_sun_vec_ECI);

    Area_m2 = SatelliteProperties.AREA_X_FACE_M2 * abs(dot(v_hat_unit,x_hat_unit_vec_ECI)) + ...
                SatelliteProperties.AREA_Y_FACE_M2 * abs(dot(v_hat_unit,y_hat_unit_vec_ECI)) + ...
                  SatelliteProperties.AREA_Z_FACE_M2 * abs(dot(v_hat_unit,z_hat_unit_vec_ECI)) + ...
                    SatelliteProperties.AREA_SOLAR_PANEL_M2 * abs(dot(v_hat_unit,solar_panel_unit_vec));

    density_kg_m3 = PhysicsConstants.RHO_0_DRAG_KG_M3 *...
                        exp(-(norm(r_ECI_meters/Units.KILOMETERS) - PhysicsConstants.R_0_DRAG_KM)/PhysicsConstants.H_DRAG_KM); % kg/m^3

    accel_drag_meters_s2 = -(1/2)*C_drag*(Area_m2/mass_kg)*density_kg_m3*vel_relative_meters_s*norm(vel_relative_meters_s); % (m2/kg)*(kg/m3)*(m/s)*(m/s) = m*m/(m*s^2) = m/s^2
end
