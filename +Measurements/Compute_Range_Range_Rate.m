function G_measurements = Compute_Range_Range_Rate(r_sat_ECI_meters, v_sat_ECI_meters_s, r_station_ECI_meters, v_station_ECI_meters_s)

    rho_vec_ECI_meters = r_sat_ECI_meters - r_station_ECI_meters; % Range vector
    rho_norm_ECI_meters = norm(rho_vec_ECI_meters); % Range magnitude

    rho_dot_ECI_meters_sec = (rho_vec_ECI_meters' * (v_sat_ECI_meters_s - v_station_ECI_meters_s)) / rho_norm_ECI_meters; % Range rate

    G_measurements = [rho_norm_ECI_meters; rho_dot_ECI_meters_sec];
end