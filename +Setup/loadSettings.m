function [S, ENV] = loadSettings(case_name, days_to_fit, scenario_config,reload_dynamics,reload_meas)

    S.case_name = case_name;
    S.days_to_fit = days_to_fit;
    S.Scenario = Setup.loadScenario(S.case_name, S.days_to_fit);

    switch scenario_config
        case 'Acceleration Testing A'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_Accel_Test();
        case 'Final - Full 3 Days'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_Final();
            S.ref_data = Tools.get_reference_data_3Days();
            S.StateCovariances = Setup.loadStateCovariances_v0();
            S.Stations = Setup.loadStations_Final();
        case 'Final - 1 Day Subset'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_Final();
            S.ref_data = Tools.get_reference_data_1Day_Subset();
            S.StateCovariances = Setup.loadStateCovariances_v0();
            S.Stations = Setup.loadStations_Final();
        case 'HW5 Subset - 6 Hours'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_HW5();
            S.ref_data = Tools.get_reference_data_HW5();
            S.StateCovariances = Setup.loadStateCovariances_v0();
            S.Stations = Setup.loadStations_HW5();
    end

    EOP_data = Tools.get_EOP_data();
    ENV.EOP_IERS = EOP_data.EOP_IERS;
    ENV.EOP_Celestrak = EOP_data.EOP_Celestrak;
    ENV.EOP_t0 = Tools.interpolate_EOP(S.IC_Sat_Epoch.epoch_date_time_UTC, ENV.EOP_IERS, ENV.EOP_Celestrak);
    ENV.time_struct = Tools.ComputeTimeSystems(S.IC_Sat_Epoch.epoch_date_time_UTC);
    ENV.EOP_IERS.MJD_days = ENV.time_struct.mjd_UTC_days;

    if reload_dynamics
        Forces.master_dynamics();
        
    end

    if reload_meas
        Measurements.master_measurements();
    end

    % TODO: Make cases for dynamic model -- switching things off and on.


end