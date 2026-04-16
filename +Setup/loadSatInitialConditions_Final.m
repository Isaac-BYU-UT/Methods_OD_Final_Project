function IC_Sat_Epoch = loadSatInitialConditions_Final()
    IC_Sat_Epoch.position_ECI_meters = [6984.45711518852;1612.2547582643;13.0925904314402] * Units.KILOMETERS; % m
    IC_Sat_Epoch.velocity_ECI_meters_per_second = [-1.67667852227336;7.26143715396544;0.259889857225218] * Units.KILOMETERS; % m/s
    IC_Sat_Epoch.epoch_date_time_UTC = datetime(2018, 3, 23, 8, 55, 3.0, 'TimeZone','UTC');
end