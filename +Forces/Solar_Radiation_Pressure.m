function accel_solar_radiation_pressure_meters_s2 = Solar_Radiation_Pressure(r_sat_rel_earth_ECI_meters, r_sun_rel_earth_ECI_meters, sat_is_illuminated)
    
    % SRP model with proper multi-face contributions.
    % Solar panel: double-gimbal, always sun-pointed
    % -Z axis: always nadir-pointed
    % X, Y faces: orthogonal to both nadir and sun directions
    
    % Unit vectors
    u_hat_sun = (r_sun_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters) / norm(r_sun_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters);
    u_hat_nadir = -r_sat_rel_earth_ECI_meters / norm(r_sat_rel_earth_ECI_meters);
    
    % Build orthogonal frame: x-axis, y-axis, z-axis (nadir)
    z_hat = u_hat_nadir;
    
    % Build x-axis perpendicular to sun and nadir
    x_axis_tmp = Tools.crossProductEquivalent(z_hat)*u_hat_sun; % cross(z_hat, u_hat_sun);
    x_hat = x_axis_tmp / norm(x_axis_tmp);
    y_hat = Tools.crossProductEquivalent(z_hat) * x_hat; %cross(z_hat, x_hat);
    y_hat = y_hat / norm(y_hat);
    
    % Face normals in spacecraft frame
    n_hat_X_pos = x_hat;
    n_hat_X_neg = -x_hat;
    n_hat_Y_pos = y_hat;
    n_hat_Y_neg = -y_hat;
    n_hat_Z_pos = -z_hat;  % +Z opposite nadir
    n_hat_Z_neg = z_hat;   % -Z is nadir-pointed
    n_hat_SP = u_hat_sun;  % solar panel always sun-pointed
    
    % Compute projected areas (dot product with sun direction)
    A_proj_X_pos = SatelliteProperties.AREA_X_FACE_M2 * max(0, n_hat_X_pos' * u_hat_sun);
    A_proj_X_neg = SatelliteProperties.AREA_X_FACE_M2 * max(0, n_hat_X_neg' * u_hat_sun);
    A_proj_Y_pos = SatelliteProperties.AREA_Y_FACE_M2 * max(0, n_hat_Y_pos' * u_hat_sun);
    A_proj_Y_neg = SatelliteProperties.AREA_Y_FACE_M2 * max(0, n_hat_Y_neg' * u_hat_sun);
    A_proj_Z_pos = SatelliteProperties.AREA_Z_FACE_M2 * max(0, n_hat_Z_pos' * u_hat_sun);
    A_proj_Z_neg = SatelliteProperties.AREA_Z_FACE_M2 * max(0, n_hat_Z_neg' * u_hat_sun);
    A_proj_SP = SatelliteProperties.AREA_SOLAR_PANEL_M2 * max(0, n_hat_SP' * u_hat_sun);
    
    % SRP pressure
    dist_sat_sun_meters = norm(r_sun_rel_earth_ECI_meters - r_sat_rel_earth_ECI_meters);
    AU_METERS = PhysicsConstants.AU_KM * Units.KILOMETERS;
    d_AU = dist_sat_sun_meters / AU_METERS;
    P_srp = PhysicsConstants.SOLAR_PRESSURE_N_M2 / d_AU^2;
    
    mass = SatelliteProperties.SATELLITE_MASS_KG;
    
    % X+ face (MLI Kapton)
    cos_phi_X_pos = max(0, n_hat_X_pos' * u_hat_sun);
    nu_X_pos = (1/3) * SatelliteProperties.MLI_KAPTON_CD;
    mu_X_pos = (1/2) * SatelliteProperties.MLI_KAPTON_CS;
    B_theta_X_pos = 2 * nu_X_pos * cos_phi_X_pos + 4 * mu_X_pos * cos_phi_X_pos^2;
    a_X_pos = -P_srp * (A_proj_X_pos / mass) * (B_theta_X_pos * n_hat_X_pos + (1 - mu_X_pos) * cos_phi_X_pos * u_hat_sun);
    
    % X- face (MLI Kapton)
    cos_phi_X_neg = max(0, n_hat_X_neg' * u_hat_sun);
    nu_X_neg = (1/3) * SatelliteProperties.MLI_KAPTON_CD;
    mu_X_neg = (1/2) * SatelliteProperties.MLI_KAPTON_CS;
    B_theta_X_neg = 2 * nu_X_neg * cos_phi_X_neg + 4 * mu_X_neg * cos_phi_X_neg^2;
    a_X_neg = -P_srp * (A_proj_X_neg / mass) * (B_theta_X_neg * n_hat_X_neg + (1 - mu_X_neg) * cos_phi_X_neg * u_hat_sun);
    
    % Y+ face (MLI Kapton)
    cos_phi_Y_pos = max(0, n_hat_Y_pos' * u_hat_sun);
    nu_Y_pos = (1/3) * SatelliteProperties.MLI_KAPTON_CD;
    mu_Y_pos = (1/2) * SatelliteProperties.MLI_KAPTON_CS;
    B_theta_Y_pos = 2 * nu_Y_pos * cos_phi_Y_pos + 4 * mu_Y_pos * cos_phi_Y_pos^2;
    a_Y_pos = -P_srp * (A_proj_Y_pos / mass) * (B_theta_Y_pos * n_hat_Y_pos + (1 - mu_Y_pos) * cos_phi_Y_pos * u_hat_sun);
    
    % Y- face (MLI Kapton)
    cos_phi_Y_neg = max(0, n_hat_Y_neg' * u_hat_sun);
    nu_Y_neg = (1/3) * SatelliteProperties.MLI_KAPTON_CD;
    mu_Y_neg = (1/2) * SatelliteProperties.MLI_KAPTON_CS;
    B_theta_Y_neg = 2 * nu_Y_neg * cos_phi_Y_neg + 4 * mu_Y_neg * cos_phi_Y_neg^2;
    a_Y_neg = -P_srp * (A_proj_Y_neg / mass) * (B_theta_Y_neg * n_hat_Y_neg + (1 - mu_Y_neg) * cos_phi_Y_neg * u_hat_sun);
    
    % Z+ face (White Paint)
    cos_phi_Z_pos = max(0, n_hat_Z_pos' * u_hat_sun);
    nu_Z_pos = (1/3) * SatelliteProperties.WHITE_PAINT_CD;
    mu_Z_pos = (1/2) * SatelliteProperties.WHITE_PAINT_CS;
    B_theta_Z_pos = 2 * nu_Z_pos * cos_phi_Z_pos + 4 * mu_Z_pos * cos_phi_Z_pos^2;
    a_Z_pos = -P_srp * (A_proj_Z_pos / mass) * (B_theta_Z_pos * n_hat_Z_pos + (1 - mu_Z_pos) * cos_phi_Z_pos * u_hat_sun);
    
    % Z- face (Germanium Kapton, nadir-pointed)
    cos_phi_Z_neg = max(0, n_hat_Z_neg' * u_hat_sun);
    nu_Z_neg = (1/3) * SatelliteProperties.GERMANIUM_KAPTON_CD;
    mu_Z_neg = (1/2) * SatelliteProperties.GERMANIUM_KAPTON_CS;
    B_theta_Z_neg = 2 * nu_Z_neg * cos_phi_Z_neg + 4 * mu_Z_neg * cos_phi_Z_neg^2;
    a_Z_neg = -P_srp * (A_proj_Z_neg / mass) * (B_theta_Z_neg * n_hat_Z_neg + (1 - mu_Z_neg) * cos_phi_Z_neg * u_hat_sun);
    
    % Solar panel (always sun-pointed, Solar Cells)
    cos_phi_SP = max(0, n_hat_SP' * u_hat_sun);
    nu_SP = (1/3) * SatelliteProperties.SOLAR_CELLS_CD;
    mu_SP = (1/2) * SatelliteProperties.SOLAR_CELLS_CS;
    B_theta_SP = 2 * nu_SP * cos_phi_SP + 4 * mu_SP * cos_phi_SP^2;
    a_SP = -P_srp * (A_proj_SP / mass) * (B_theta_SP * n_hat_SP + (1 - mu_SP) * cos_phi_SP * u_hat_sun);
    
    % Total acceleration
    accel_total = a_X_pos + a_X_neg + a_Y_pos + a_Y_neg + a_Z_pos + a_Z_neg + a_SP;
    accel_solar_radiation_pressure_meters_s2 = accel_total * sat_is_illuminated;

end