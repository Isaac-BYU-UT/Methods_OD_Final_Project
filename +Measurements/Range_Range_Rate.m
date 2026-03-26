function G = range_range_rate_sym(r_sat_ECI_km, v_sat_ECI_km_s, r_station_ECI_km, v_station_ECI_km_s)

    rho_vec_ECI_km = r_sat_ECI_km - r_station_ECI_km; % Range vector
    rho_norm_ECI_km = norm(rho_vec_ECI_km); % Range magnitude

    rho_dot_ECI_km_sec = (rho_vec_ECI_km' * (v_sat_ECI_km_s - v_station_ECI_km_s)) / rho_norm_ECI_km; % Range rate

    G = [rho_norm_ECI_km; rho_dot_ECI_km_sec];
end