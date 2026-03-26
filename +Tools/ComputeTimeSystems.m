function Time_Struct = ComputeTimeSystems(utc_datetime, delta_AT_sec, delta_UT1_sec)
    % --- CONSTANTS (Explicit Units) ---
    JD_J2000_days    = 2451545.0;   % J2000.0 epoch
    MJD_OFFSET_days  = 2400000.5;   % Offset for Modified Julian Date
    SEC_PER_DAY      = 86400;       
    DAYS_PER_CENTURY = 36525;       

    % --- 1. UTC & MJD (Always Calculated) ---
    Time_Struct.jd_UTC_days  = juliandate(utc_datetime);
    Time_Struct.mjd_UTC_days = Time_Struct.jd_UTC_days - MJD_OFFSET_days;
    
    % --- 2. UT1 (Universal Time 1) ---
    % Only calculate if delta_UT1_sec is provided (at least 3 arguments)
    if nargin >= 3 && ~isempty(delta_UT1_sec)
        Time_Struct.jd_UT1_days = Time_Struct.jd_UTC_days + (delta_UT1_sec / SEC_PER_DAY);
        
        % Calculation for Midnight (0h) UT1
        Time_Struct.jd_UT1_0h_days = floor(Time_Struct.jd_UT1_days + 0.5) - 0.5;
        
        % Elapsed seconds since midnight UT1
        Time_Struct.sec_past_midnight_UT1 = (Time_Struct.jd_UT1_days - Time_Struct.jd_UT1_0h_days) * SEC_PER_DAY;
        
        % Julian Centuries at 0h UT1
        Time_Struct.t_UT1_0h_centuries = (Time_Struct.jd_UT1_0h_days - JD_J2000_days) / DAYS_PER_CENTURY;
    end

    % --- 3. TT (Terrestrial Time) ---
    % Only calculate if delta_AT_sec is provided (at least 2 arguments)
    if nargin >= 2 && ~isempty(delta_AT_sec)
        % TT = UTC + delta_AT (TAI-UTC) + 32.184s (TT-TAI)
        TT_offset_from_UTC_sec = delta_AT_sec + 32.184;
        Time_Struct.jd_TT_days = Time_Struct.jd_UTC_days + (TT_offset_from_UTC_sec / SEC_PER_DAY);
        
        % Julian Centuries from J2000.0
        Time_Struct.t_TT_centuries = (Time_Struct.jd_TT_days - JD_J2000_days) / DAYS_PER_CENTURY;
    end

end