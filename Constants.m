classdef Constants
    properties (Constant)
        KM_TO_METERS      = 1000;
        METERS_TO_KM      = 1/1000;
        MU_EARTH_KM3_S2   = 398600.4415
        R_EARTH_KM        = 6378.1363
        MU_SUN_KM3_S2     = 132712440018
        AU_KM             = 149597870.7
        MU_MOON_KM3_S2    = 4902.800066
        ECC_EARTH         = 0.081819221456
        OMEGA_EARTH_RAD_S = 7.292115146706979e-5
        RHO_0_DRAG_KG_M3  = 3.614e-13
        R_0_DRAG_KM        = (700000.0 * (1/1000) + 6378.1363) % R_EARTH_KM + 700000 m * METERS_TO_KM
        H_DRAG_KM          = 88667.0 * (1/1000) % We want this in KM!
        
        J2_EARTH          =  1.08262617385222e-03 % from Vallado D-1 (EGM-08)
        J3_EARTH          = -2.53241051856772e-06 % from Vallado D-1 (EGM-08)
        J4_EARTH          = -1.61989759991697e-06 % from Vallado D-1 (EGM-08)

        ARCSEC_TO_RAD     = pi/648000 % Number of arseconds per radian
        ARCSEC_TO_DEG     = 1/3600 % Number of arcseconds per degree
        MILLI_TO_NOM      = 1/1000 % Number of nominal unit per MILLI
        SEC_IN_SOLAR_DAY  = 24 * 60 * 60 %24 hrs * 60 min * 60 sec
        
        SATTELITE_MASS_KG   = 2000 % Mass of satellite for drag calculations, in kg
        % TODO -- this will be replaced later
        SATELLITE_AREA_M2   = 20 % Cross-sectional area of satellite for drag calculations, in m^2

        % For validation purposes only
        % R_EARTH_KM        = 6378.145 % From HW2
        % MU_EARTH_KM3_S2   = 398600.4 % From HW2
    end
end

%% TODO: DOUBLE CHECK DECIMAL PRECISION ON THESE!!