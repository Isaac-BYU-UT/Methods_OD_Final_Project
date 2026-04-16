function IC_Sat_Epoch = loadSatInitialConditions_HW5()
    IC_Sat_Epoch.position_ECI_meters = [6990077.798814194;1617465.311978378;22679.810569245355]; % m
    IC_Sat_Epoch.velocity_ECI_meters_per_second = [-1675.13972506056;7273.72441330686;252.688512916741]; % m/s
    IC_Sat_Epoch.epoch_date_time_UTC = datetime(2018, 2, 1, 5, 0, 0.0, 'TimeZone','UTC');
end