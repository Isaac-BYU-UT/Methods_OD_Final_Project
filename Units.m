classdef Units
    properties (Constant)
        KILOMETERS = 1000;
        HOURS = 3600;
        MILIMETERS = 1/1000;

        ARCSEC_TO_RAD     = pi/648000 % Number of arseconds per radian
        ARCSEC_TO_DEG     = 1/3600 % Number of arcseconds per degree
        RAD_TO_DEG        = 180/pi % Number of degrees per radian
        DEG_TO_RAD        = pi/180 % Number of radians per degree.
        MILLI_TO_NOM      = 1/1000 % Number of nominal unit per MILLI
        SEC_IN_SOLAR_DAY  = 24 * 60 * 60 %24 hrs * 60 min * 60 sec
        
    end
end