function IC = initial_conditions_Phase1_EXCLUDE_RANGE_RATE()

IC.name = 'Final Phase 1 Exclude Range Rate';
IC.description = 'Final Phase 1 Exclude Range Rate';    

IC.plotting_histograms = false;

%  EPOCH
% ========================
IC.UTC_epoch = datetime(2018, 3, 23, 8, 55, 3.0, 'TimeZone','UTC');
IC.time_struct = Tools.ComputeTimeSystems(IC.UTC_epoch);

%  SATELLITE INITIAL STATE
% ========================
IC.r_sat_ECI_km = [6984.45711518852;1612.2547582643;13.0925904314402] ...
                   ; % position given in km!

IC.v_sat_ECI_km_s = [-1.67667852227336;7.26143715396544;0.259889857225218] ...
                    ; % veloctiy given in km/sec!

IC.C_drag = 1.88;

% State Covariance
position_std_km = 10; % km standard deviation
velocity_std_km_s = 10 * Constants.METERS_TO_KM; % 10 m/s to km/s
drag_coefficient_std = .05; % THIS IS JUST A GUESS!
IC.P_Covariance_States = diag([position_std_km^2, position_std_km^2, position_std_km^2, ...
                            velocity_std_km_s^2, velocity_std_km_s^2, velocity_std_km_s^2, ...
                            drag_coefficient_std^2]);
IC.Sigma_Accel_km_s2 = 10e-9; % TODO: TUNE THIS GUESS! % If we go larger than 10e-8, we get problems!

%  STATIONS (ECEF)
% ========================
station_names = {'Atoll', 'Diego Garcia', 'Arecibo'};

raw_coords_ECEF_m = [
    -6143584,  1364250,  1033743;
     1907295,  6030810,  -817119;
     2390310, -5564341,  1994578
];

IC.Stations = struct();

% Each station has zero mean noise, range_std = 5 m and range-rate std = 1.0 mm/sec
% TODO : For future version, each station will have a different covariance.

Almost_Infinite = 10e18;

Range_STD_Atoll_km = 10 * Constants.METERS_TO_KM;
Range_Rate_STD_Atoll_km_s = 0.5 * Constants.METERS_TO_KM / 1000 * Almost_Infinite; % mm/s to km/s

Range_STD_Diego_Garcia_km = 5 * Constants.METERS_TO_KM;
Range_Rate_STD_Diego_Garcia_km_s = 1.0 * Constants.METERS_TO_KM / 1000 * Almost_Infinite; % mm/s to km/s

Range_STD_Arecibo_km = 10 * Constants.METERS_TO_KM;
Range_Rate_STD_Arecibo_km_s = 0.5 * Constants.METERS_TO_KM / 1000 * Almost_Infinite; % mm/s to km/s

for i = 1:length(station_names)
    IC.Stations(i).Name = station_names{i};
    IC.Stations(i).r_ECEF_km = raw_coords_ECEF_m(i, :)' * Constants.METERS_TO_KM;
end

% Add covariance matrix for each station
IC.Stations(1).Covariance = diag([Range_STD_Atoll_km^2, Range_Rate_STD_Atoll_km_s^2]);
IC.Stations(2).Covariance = diag([Range_STD_Diego_Garcia_km^2, Range_Rate_STD_Diego_Garcia_km_s^2]);
IC.Stations(3).Covariance = diag([Range_STD_Arecibo_km^2, Range_Rate_STD_Arecibo_km_s^2]);

%  TIME SETTINGS (ODE)
% ========================
IC.time.end = 21600; % sec
IC.time.dt  = 60;    % sec
IC.time.tspan = 0:IC.time.dt:IC.time.end;

end