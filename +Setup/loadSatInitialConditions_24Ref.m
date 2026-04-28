function IC_Sat_Epoch = loadSatInitialConditions_24Ref()
    
    State = load("ref_data/24hourStates.mat").State;
    
    r_ref_init_km = State(1,1:3);
    v_ref_init_km = State(1,4:6);
    
    r_init_m = r_ref_init_km(:) * Units.KILOMETERS;
    v_init_m_s = v_ref_init_km(:) * Units.KILOMETERS;

    IC_Sat_Epoch.position_ECI_meters = r_init_m; % m
    IC_Sat_Epoch.velocity_ECI_meters_per_second = v_init_m_s; % m/s
    IC_Sat_Epoch.epoch_date_time_UTC = datetime(2018, 3, 23, 8, 55, 3.0, 'TimeZone','UTC');
end