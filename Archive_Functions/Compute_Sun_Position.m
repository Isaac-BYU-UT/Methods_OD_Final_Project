function r_vec_Sun_MOD_km = Compute_Sun_Position(T_UT1) % T_UT1 is (JD_UT1 - J2000.0)/(Days per JC)
    
    % SEE VALLADO 5.1

    % Textbook constants (updated for precision)
    mean_longitude_sun_DEG = 280.460 + 36000.771 * T_UT1;
    mean_anomaly_sun_DEG = 357.52910922 + 35999.05034 * T_UT1;
    
    % Normalization
    mean_longitude_sun_DEG = mod(mean_longitude_sun_DEG, 360);
    mean_anomaly_sun_DEG = mod(mean_anomaly_sun_DEG, 360);
    
    % Ecliptic Longitude
    ecliptic_long_deg = mean_longitude_sun_DEG + ...
                        1.914666471 * sind(mean_anomaly_sun_DEG) + ...
                        0.019994643 * sind(2 * mean_anomaly_sun_DEG);
               
    % Distance in AU
    r_sun_AU = 1.000140612 - 0.016708617 * cosd(mean_anomaly_sun_DEG) - ...
               0.000139589 * cosd(2 * mean_anomaly_sun_DEG);
          
    % Obliquity of the Ecliptic
    eps_deg = 23.439291 - 0.01300427 * T_UT1;
    
    % Form Vector (MOD Frame)
    r_vector_Sun_MOD_AU = r_sun_AU .* [ ...
                            cosd(ecliptic_long_deg); ...
                            cosd(eps_deg) * sind(ecliptic_long_deg); ...
                            sind(eps_deg) * sind(ecliptic_long_deg) ...
                         ];

    r_vec_Sun_MOD_km = r_vector_Sun_MOD_AU * Constants.AU_KM;

end