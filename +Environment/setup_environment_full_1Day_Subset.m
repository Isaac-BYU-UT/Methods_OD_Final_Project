function ENV = setup_environment_full_1Day_Subset(IC)

%% Load Data
ENV.ref_data = Tools.get_reference_data_1Day_Subset();
EOP_data = Tools.get_EOP_data();

ENV.EOP_IERS = EOP_data.table_EOP_IERS;
ENV.EOP_Celestrak = EOP_data.table_EOP_Celestrak;

%% Functions
[ENV.AccelFcn, ENV.AmatrixFcn] = Forces.master_dynamics();
[ENV.MeasFcn, ENV.HmatrixFcn] = Measurements.master_measurements();

%% EOP at t0
ENV.EOP_t0 = Tools.interpolate_EOP(IC.UTC_epoch, "IERS");

%% Convert Stations to ECI at t0
for i = 1:length(IC.Stations)

    [r_ECI, v_ECI] = Tools.ECEF_ECI_Converter( ...
        IC.Stations(i).r_ECEF_km, ...
        zeros(3,1), ...
        IC.UTC_epoch, ...
        'ECEF_to_ECI', ...
        ENV.EOP_t0);

    ENV.Stations(i).Name = IC.Stations(i).Name;
    ENV.Stations(i).r_ECEF_km = IC.Stations(i).r_ECEF_km;
    ENV.Stations(i).r_ECI_km = r_ECI;
    ENV.Stations(i).v_ECI_km_s = v_ECI;
end

%% Sun / Moon
[ENV.r_sun_ECI_km, ~] = Forces.Vallado_sunPositionLowPrecision(IC.time_struct.jd_UTC_days);
[ENV.r_moon_ECI_km, ~] = Forces.Vallado_moonPositionLowPrecision(IC.time_struct.jd_UTC_days);

ENV.sunLOS = Forces.Vallado_sunLOS(IC.r_sat_ECI_km, ENV.r_sun_ECI_km);

end