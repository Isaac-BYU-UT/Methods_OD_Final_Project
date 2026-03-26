function EOP_Data = get_EOP_data()

    persistent EOP_Data_cached

    if isempty(EOP_Data_cached)
        fprintf('Loading Earth Orientation Parameters (EOP) data from file...\n');

        %% Import EOP from IERS Bulletin A and B (finals.all.csv)
        % Column names: MJD;Year;Month;Day;Type;x_pole;sigma_x_pole;y_pole;sigma_y_pole;x_rate;sigma_x_rate;y_rate;sigma_y_rate;Type;UT1-UTC;sigma_UT1-UTC;LOD;sigma_LOD;Type;dPsi;sigma_dPsi;dEpsilon;sigma_dEpsilon;dX;sigma_dX;dY;sigma_dY;Type;bulB/x_pole;bulB/y_pole;Type;bulB/UT-UTC;Type;bulB/dPsi;bulB/dEpsilon;bulB/dX;bulB/dY
        filename = 'EOP_data/finals.all.csv';
        opts = detectImportOptions(filename, 'Delimiter', ';');
        opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', '');     % Tell MATLAB that empty fields should be treated as missing
        opts.VariableNamingRule = 'preserve';
        table_EOP_IERS = readtable(filename, opts);

        % Rename variables with clear names + units
        table_EOP_IERS.Properties.VariableNames = { ...
            'MJD_days', ...
            'Year', ...
            'Month', ...
            'Day', ...
            'PM_flag_A', ...
            'x_pole_arcsec', ...
            'sigma_x_pole_arcsec', ...
            'y_pole_arcsec', ...
            'sigma_y_pole_arcsec', ...
            'x_rate_arcsec_per_day', ...
            'sigma_x_rate_arcsec_per_day', ...
            'y_rate_arcsec_per_day', ...
            'sigma_y_rate_arcsec_per_day', ...
            'UT1_flag_A', ...
            'UT1_minus_UTC_sec', ...
            'sigma_UT1_minus_UTC_sec', ...
            'LOD_millisec', ...
            'sigma_LOD_millisec', ...
            'Nutation_flag_A', ...
            'dPsi_milli_arcsec', ...
            'sigma_dPsi_milli_arcsec', ...
            'dEpsilon_milli_arcsec', ...
            'sigma_dEpsilon_milli_arcsec', ...
            'dX_milli_arcsec', ...
            'sigma_dX_milli_arcsec', ...
            'dY_milli_arcsec', ...
            'sigma_dY_milli_arcsec', ...
            'PM_flag_B', ...
            'x_pole_B_arcsec', ...
            'y_pole_B_arcsec', ...
            'UT1_flag_B', ...
            'UT1_minus_UTC_B_sec', ...
            'Nutation_flag_B', ...
            'dPsi_B_milli_arcsec', ...
            'dEpsilon_B_milli_arcsec', ...
            'dX_B_milli_arcsec', ...
            'dY_B_milli_arcsec' ...
        };

        table_EOP_IERS.Date = datetime( ...
                                        table_EOP_IERS.Year, ...
                                        table_EOP_IERS.Month, ...
                                        table_EOP_IERS.Day, ...
                                        'TimeZone', 'UTC' );

        table_EOP_IERS = movevars(table_EOP_IERS, 'Date', 'Before', 'MJD_days');

        % OP IAU1980 from https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
        % Year, Month, Day, Modified Julian Date

        % PM-x [arcsec], error_PM-x [arcsec],
        % PM-y [arcsec], error_PM-y [arcsec],
        % UT1-UTC [seconds], error_UT1-UTC [seconds],
        % LOD [milliseconds], error_LOD [milliseconds],
        % dPsi [milliarcsec], error_dPsi [milliarcsec],
        % dEps [milliarcsec], error_dEps [milliarcsec]
        % (from Bulletin A)

        % PM-x [arcsec], PM-y [arcsec], UT1-UTC [seconds],
        % dPsi [milliarcsec], dEps [milliarcsec]
        % (from Bulletin B)

        %         Col.#    Format  Quantity
        % -------  ------  -------------------------------------------------------------
        % 1-2      I2      year (to get true calendar year, add 1900 for MJD<=51543 or add 2000 for MJD>=51544)
        % 3-4      I2      month number
        % 5-6      I2      day of month
        % 7        X       [blank]
        % 8-15     F8.2    fractional Modified Julian Date (MJD UTC)
        % 16       X       [blank]
        % 17       A1      IERS (I) or Prediction (P) flag for Bull. A polar motion values
        % 18       X       [blank]
        % 19-27    F9.6    Bull. A PM-x (sec. of arc)
        % 28-36    F9.6    error in PM-x (sec. of arc)
        % 37       X       [blank]
        % 38-46    F9.6    Bull. A PM-y (sec. of arc)
        % 47-55    F9.6    error in PM-y (sec. of arc)
        % 56-57    2X      [blanks]
        % 58       A1      IERS (I) or Prediction (P) flag for Bull. A UT1-UTC values
        % 59-68    F10.7   Bull. A UT1-UTC (sec. of time)
        % 69-78    F10.7   error in UT1-UTC (sec. of time)
        % 79       X       [blank]
        % 80-86    F7.4    Bull. A LOD (msec. of time) -- NOT ALWAYS FILLED
        % 87-93    F7.4    error in LOD (msec. of time) -- NOT ALWAYS FILLED
        % 94-95    2X      [blanks]
        % 96       A1      IERS (I) or Prediction (P) flag for Bull. A nutation values
        % 97       X       [blank]
        % 98-106   F9.3    Bull. A dPSI (msec. of arc)
        % 107-115  F9.3    error in dPSI (msec. of arc)
        % 116      X       [blank]
        % 117-125  F9.3    Bull. A dEPSILON (msec. of arc)
        % 126-134  F9.3    error in dEPSILON (msec. of arc)
        % 135-144  F10.6   Bull. B PM-x (sec. of arc)
        % 145-154  F10.6   Bull. B PM-y (sec. of arc)
        % 155-165  F11.7   Bull. B UT1-UTC (sec. of time)
        % 166-175  F10.3   Bull. B dPSI (msec. of arc)
        % 176-185  F10.3   Bull. B dEPSILON (msec. of arc)

        %% Import Data from Celestrak EOP-All.txt

        % Celestrak Dataset for EOP: https://celestrak.org/SpaceData/EOP-All.txt
        % 001-004	Year
        % 006-007	Month (01-12)
        % 009-010	Day
        % 012-016	Modified Julian Date (Julian Date at 0h UTC minus 2400000.5)
        % 018-026	x (arc seconds)
        % 028-036	y (arc seconds)
        % 038-047	UT1-UTC (seconds)
        % 049-058	Length of Day (seconds)
        % 060-068	δΔψ (arc seconds)
        % 070-078	δΔε (arc seconds)
        % 080-088	δX (arc seconds)
        % 090-098	δY (arc seconds)
        % 100-102	Delta Atomic Time, TAI-UTC (seconds)

        filename = 'EOP_data/EOP-All.txt';

        % Open file
        fid = fopen(filename, 'r');

        % Move to "BEGIN OBSERVED"
        while true
            line = fgetl(fid);
            if contains(line, 'BEGIN OBSERVED')
                break;
            end
        end

        % Read numeric data until "END PREDICTED"
        data = [];
        while true
            line = fgetl(fid);
            
            if ~ischar(line) || contains(line, 'END PREDICTED')
                break;
            end
            
            % Parse numeric line
            nums = sscanf(line, '%f');
            
            % Only keep valid rows (should have 13 values)
            if numel(nums) == 13
                data(end+1, :) = nums'; %#ok<SAGROW>
            end
        end

        fclose(fid);

        % Convert to table with labeled columns
        table_EOP_Celestrak = array2table(data, 'VariableNames', { ...
            'Year', ...
            'Month', ...
            'Day', ...
            'MJD_days', ...
            'x_pole_arcsec', ...
            'y_pole_arcsec', ...
            'UT1_minus_UTC_sec', ...
            'LOD_sec', ...
            'dPsi_arcsec', ...
            'dEpsilon_arcsec', ...
            'dX_arcsec', ...
            'dY_arcsec', ...
            'TAI_minus_UTC_sec' ...
        });

        % Add Extra Columns to Match IERS Data + Convert Units
        table_EOP_Celestrak.Date = datetime(table_EOP_Celestrak.Year, table_EOP_Celestrak.Month, table_EOP_Celestrak.Day,'TimeZone','UTC');
        table_EOP_Celestrak.LOD_millisec = table_EOP_Celestrak.LOD_sec * 1000; % Convert LOD to milliseconds for consistency with IERS data
        table_EOP_Celestrak.dPsi_milli_arcsec = table_EOP_Celestrak.dPsi_arcsec * 1000; % Convert dPsi to milli-arcseconds for consistency with IERS data
        table_EOP_Celestrak.dEpsilon_milli_arcsec = table_EOP_Celestrak.dEpsilon_arcsec * 1000; % Convert dEpsilon to milli-arcseconds for consistency with IERS data
        table_EOP_Celestrak.dX_milli_arcsec = table_EOP_Celestrak.dX_arcsec * 1000; % Convert dX to milli-arcseconds for consistency with IERS data
        table_EOP_Celestrak.dY_milli_arcsec = table_EOP_Celestrak.dY_arcsec * 1000; % Convert dY to milli-arcseconds for consistency with IERS data

        % Move Date to front
        table_EOP_Celestrak = movevars(table_EOP_Celestrak, 'Date', 'Before', 'Year');

        %% Import Data from nut80.dat
        % Okay Now, I need to use the Dat80.m to get my nut80 nutation parameters.
        nut80 = readtable("EOP_data/nut80.dat");

        % Table nut80.dat has sample coefficients for integers (an1i− an5i) and the real coefficients (Ai, Bi, Ci, Di).
        nut80.Properties.VariableNames = {'a_n1','a_n2','a_n3','a_n4','a_n5','A_i_1e_4_arcsec','B_i_1e_4_arcsec','C_i_1e_4_arcsec','D_i_1e_4_arcsec','i_table'};

        % TABLE D-6. Largest 1980 IAU Theory of Nutation Coefficients, Epoch J2000.  The units for
        % the longitude terms (Ai, Bi) and the obliquity terms (Ci, Di) are 0.0001'' per Julian
        % century. Data is from Seidelmann (1992:112–113) and McCarthy (1996). There are    
        % 106 terms in the complete theory. 

        unit_conversion = .0001 * Constants.ARCSEC_TO_DEG; % .0001" per Julian Century * 

        nut80.A_i_deg = nut80.A_i_1e_4_arcsec * unit_conversion;
        nut80.B_i_deg = nut80.B_i_1e_4_arcsec * unit_conversion;
        nut80.C_i_deg = nut80.C_i_1e_4_arcsec * unit_conversion;
        nut80.D_i_deg = nut80.D_i_1e_4_arcsec * unit_conversion;
     
        %% Set Equal to Cashed Data
        EOP_Data_cached = struct('table_EOP_IERS', table_EOP_IERS, 'table_EOP_Celestrak', table_EOP_Celestrak, 'nut80', nut80);

    end
        EOP_Data = EOP_Data_cached;

    
    
end