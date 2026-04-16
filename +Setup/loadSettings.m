function [S, ENV] = loadSettings(case_name, scenario_config, reload_dynamics, reload_meas)
    
if nargin < 3
    reload_dynamics = false;
end
if nargin < 4
    reload_meas = false;
end

S.case_name = case_name;
S.scenario_config = scenario_config;
%% Load in Initial Conditions
    switch scenario_config
        case 'ECI_ECEF_test'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_ECI_ECEF_Test();
        case 'Accel' %'Acceleration Testing A'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_Accel_Test();
        case 'Final_3D' %'Final - Full 3 Days'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_Final();
            S.ref_data = Tools.get_reference_data_3Days();
            S.StateCovariances = Setup.loadStateCovariances_v0();
            S.Stations = Setup.loadStations_Final();
        case 'Final_1D' %'Final - 1 Day Subset'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_Final();
            S.ref_data = Tools.get_reference_data_1Day_Subset();
            S.StateCovariances = Setup.loadStateCovariances_v0();
            S.Stations = Setup.loadStations_Final();
        case 'HW5' %'HW5 Subset - 6 Hours'
            S.IC_Sat_Epoch = Setup.loadSatInitialConditions_HW5();
            S.ref_data = Tools.get_reference_data_HW5();
            S.StateCovariances = Setup.loadStateCovariances_v0();
            S.Stations = Setup.loadStations_HW5();
    end
%% Load in EOP Stuff
    EOP_data = Tools.get_EOP_data();
    ENV.EOP_IERS = EOP_data.EOP_IERS;
    ENV.EOP_IERS.t = ENV.EOP_IERS.MJD_days;
    ENV.EOP_IERS.Y = [ENV.EOP_IERS.x_pole_arcsec, ENV.EOP_IERS.y_pole_arcsec, ENV.EOP_IERS.dPsi_milli_arcsec, ...
                  ENV.EOP_IERS.dEpsilon_milli_arcsec, ENV.EOP_IERS.UT1_minus_UTC_sec, ENV.EOP_IERS.LOD_millisec];
    ENV.EOP_Celestrak = EOP_data.EOP_Celestrak;

    % Trim EOP data to just the year of our epoch for efficiency (and to avoid interpolation issues at the boundaries)
    EOP_year_mask = year(S.IC_Sat_Epoch.epoch_date_time_UTC) == ENV.EOP_IERS.Year;
    ENV.EOP_IERS.t = ENV.EOP_IERS.t(EOP_year_mask);
    ENV.EOP_IERS.Y = ENV.EOP_IERS.Y(EOP_year_mask,:);

    % Trim Celestrak data to the same year for consistency
    Celestrak_year_mask = year(S.IC_Sat_Epoch.epoch_date_time_UTC) == ENV.EOP_Celestrak.Year;
    ENV.EOP_Celestrak.MJD_days = ENV.EOP_Celestrak.MJD_days(Celestrak_year_mask);
    ENV.EOP_Celestrak.TAI_minus_UTC_sec = ENV.EOP_Celestrak.TAI_minus_UTC_sec(Celestrak_year_mask);

    ENV.EOP_t0 = Tools.interpolate_EOP(S.IC_Sat_Epoch.epoch_date_time_UTC, ENV.EOP_IERS, ENV.EOP_Celestrak);
    ENV.time_struct_epoch = Tools.ComputeTimeSystems(S.IC_Sat_Epoch.epoch_date_time_UTC);


%% Load Dyanmics and Measurement Functions
    if reload_dynamics
        Forces.master_dynamics();
    end

    if reload_meas
        Measurements.master_measurements();
    end
%%  Load Scenario Params
    S.Scenario.range_on = true;
    S.Scenario.range_rate_on = true;
    S.Scenario.Atoll_on = true;
    S.Scenario.Diego_Garcia_on = true;
    S.Scenario.Arecibo_on = true;
    
    switch case_name
        case 'A' % Range Only, All Sensors
            S.Scenario.range_rate_on = false;
            return;
    
        case 'B' % Range-Rate Only, All Sensors
            S.Scenario.range_on = false;
            return;
    
        case 'C' % Atoll Only, All Data Types
            S.Scenario.Diego_Garcia_on = false;
            S.Scenario.Arecibo_on = false;
            return;
    
        case 'D' % Diego Garcia Only, All Data Types
            S.Scenario.Atoll_on = false;
            S.Scenario.Arecibo_on = false;
            return;
    
        case 'E' % Arecibo Only, All Data Types
            S.Scenario.Atoll_on = false;
            S.Scenario.Diego_Garcia_on = false;
            return;
    
        case 'F' % Fit the Long Arc, All Data, All Stations
            return;
    
         % TODO: FIGURE THIS ONE OUT!
        case 'G' % Fit the Short Arc, Only The Last Day of Data for All Sensors
            % TODO: Need to implement the "short arc" part of this, but for now just return the same as F
            return;

    end


end