classdef SatelliteProperties
    properties (Constant)
        SATELLITE_MASS_KG = 2000; % Mass of satellite for drag calculations, in kg
       
        % --- Spacecraft Geometry and Optical Properties ---
        % Areas in m^2
        AREA_X_FACE_M2 = 6.0;
        AREA_Y_FACE_M2 = 8.0;
        AREA_Z_FACE_M2 = 12.0; % Same for +Z and -Z
        AREA_SOLAR_PANEL_M2 = 15.0;

        % Optical Coefficients (Cd = Diffuse, Cs = Specular)
        % MLI Kapton (+X, -X, +Y, -Y)
        MLI_KAPTON_CD = 0.04;
        MLI_KAPTON_CS = 0.59;

        % White Paint (+Z)
        WHITE_PAINT_CD = 0.80;
        WHITE_PAINT_CS = 0.04;

        % Germanium Kapton (-Z)
        GERMANIUM_KAPTON_CD = 0.28;
        GERMANIUM_KAPTON_CS = 0.18;

        % Solar Cells
        SOLAR_CELLS_CD = 0.04;
        SOLAR_CELLS_CS = 0.04;

        % PLACEHOLDER FOR NOW -- ASSUME THAT C_DRAG IS ALWAYS CONSTANT!!
        C_Drag = 1.88;
    end
end