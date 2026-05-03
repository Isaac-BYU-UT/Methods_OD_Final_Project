function IC_Sat_Epoch = loadSatInitialConditions_Final_Batched()
    IC_Sat_Epoch.position_ECI_meters = [ 6.978647134050959e+06; 1.616545991696636e+06; 1.949657928133458e+04]; % m
    IC_Sat_Epoch.velocity_ECI_meters_per_second = [-1.662986917069318e+03; 7.260785572218775e+03; 2.706036897075662e+02]; % m/s
    IC_Sat_Epoch.epoch_date_time_UTC = datetime(2018, 3, 23, 8, 55, 3.0, 'TimeZone','UTC');
end