function IC_Sat_Epoch = loadSatInitialConditions_ECI_ECEF_Test()
    IC_Sat_Epoch.position_ECI_meters = []; % m

    IC_Sat_Epoch.velocity_ECI_meters_per_second = []; % m/s

    IC_Sat_Epoch.epoch_date_time_UTC = datetime(2004, 4, 6, 7, 51, 28.386009, 'TimeZone','UTC');
end