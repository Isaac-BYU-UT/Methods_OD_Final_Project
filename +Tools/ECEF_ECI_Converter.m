function [r_converted_km, v_converted_km_sec] = ECEF_ECI_Converter(r_original_km, v_original_km_sec, UTC_date_time, conversion_type, EOP_params)

time_struct = Tools.ComputeTimeSystems(UTC_date_time, EOP_params.delta_AT_sec, EOP_params.UT1_UTC_sec);
arcsec2rad = pi/648000; % Number of arseconds per radian
arcsec2deg = 1/3600; % Number of arcseconds per degree 
t_TT = time_struct.t_TT_centuries;
%% P Procession Matrix
zeta_rad  = (2306.2181*t_TT + 0.30188*t_TT^2 + 0.017998*t_TT^3) * Constants.ARCSEC_TO_RAD;
theta_rad = (2004.3109*t_TT - 0.42665*t_TT^2 - 0.041833*t_TT^3) * Constants.ARCSEC_TO_RAD;
z_rad     = (2306.2181*t_TT + 1.09468*t_TT^2 + 0.018203*t_TT^3) * Constants.ARCSEC_TO_RAD;

P = rot3_rad(zeta_rad)*rot2_rad(-theta_rad)*rot3_rad(z_rad);
%% N Nutation Matrix

% Pull out nut80 params
ECI_ECEF_Transform_Data = Tools.get_EOP_data(); % Includes table_EOP_IERS, table_EOP_Celestrak, and nut80
nut80 = ECI_ECEF_Transform_Data.nut80;

% Using Eq 3-68
epsilon_mean_1980_arcsec = 84381.448 - 46.8150*t_TT - 0.00059*t_TT^2 + 0.001813*t_TT^3; % in arcseconds
epsilon_mean_1980_deg = 23.439291 - 0.0130042*t_TT - (1.64e-7)*t_TT^2 + (5.04e-7)*t_TT^3; % in deg

% Delaunay Parameters (FROM ERRATA!)
r_deg = 360; % 360 degrees

% Mean anomaly of the moon
M_moon_deg = 134.96298139 ... 
    + (1325*r_deg + 198.8673981)*t_TT ...
    + 0.0086972*t_TT^2 ...
    + 1.78e-5*t_TT^3;

% Mean anomaly of the sun
M_sun_deg = 357.52772333 ... 
    + (99*r_deg + 359.0503400)*t_TT ...
    - 0.0001603*t_TT^2 ...
    - 3.3e-6*t_TT^3;

% Mean Argument of Latitude of the Moon
u_M_moon_deg = 93.27191028 ... 
    + (1342*r_deg + 82.0175381)*t_TT ...
    - 0.0036825*t_TT^2 ...
    + 3.1e-6*t_TT^3;

% Mean Elongation of the Moon from the Sun
D_deg = 297.85036306 ...    
    + (1236*r_deg + 307.1114800)*t_TT ...
    - 0.0019142*t_TT^2 ...
    + 5.3e-6*t_TT^3;

% Longitude of the Ascending Node of the Moon
Omega_deg = 125.04452222 ... 
    - (5*r_deg + 134.1362608)*t_TT ...
    + 0.0020708*t_TT^2 ...
    + 2.2e-6*t_TT^3;

M_moon_deg = mod(M_moon_deg,360);
M_sun_deg  = mod(M_sun_deg,360);
u_M_moon_deg = mod(u_M_moon_deg,360);
D_deg      = mod(D_deg,360);
Omega_deg  = mod(Omega_deg,360);

num_nutation_params = height(nut80);

% Initialize Values
DELTA_psi_1980_deg = 0;
DELTA_epsilon_1980_deg = 0;

for i = 1:num_nutation_params % USing Vallado's Eq. 3-83
    a_p_i_deg = nut80.a_n1(i)*M_moon_deg + nut80.a_n2(i)*M_sun_deg + nut80.a_n3(i)*u_M_moon_deg + nut80.a_n4(i)*D_deg + nut80.a_n5(i)*Omega_deg;
    a_p_i_deg = mod(a_p_i_deg,360);
    DELTA_psi_1980_deg = DELTA_psi_1980_deg + (nut80.A_i_deg(i) + nut80.B_i_deg(i)*t_TT)*sind(a_p_i_deg);         % Note the sind, since all of our vars are in degrees
    DELTA_epsilon_1980_deg = DELTA_epsilon_1980_deg + (nut80.C_i_deg(i) + nut80.D_i_deg(i)*t_TT)*cosd(a_p_i_deg); % Note the cosd, since all of our vars are in degrees
end

% Eq 3-84
delta_psi_final_deg = DELTA_psi_1980_deg + EOP_params.small_d_delta_psi_1980_deg;
delta_epsilon_final_deg = DELTA_epsilon_1980_deg + EOP_params.small_d_delta_epsilon_1980_deg;

% Eq 3-85
epsilon_1980_final_deg = epsilon_mean_1980_deg + delta_epsilon_final_deg; % Dont use the small correction factor here

N = rot1_deg(-epsilon_mean_1980_deg)*rot3_deg(delta_psi_final_deg)*rot1_deg(epsilon_1980_final_deg);
%% S Siderial Time Matrix 

t_GMST = time_struct.t_UT1_0h_centuries;

THETA_GMST_1982_deg_0h = ... % Vallado Eq 3-45 degree version
      100.4606184 ...
    + 36000.77005361 * t_GMST ...
    + 0.00038793 * t_GMST^2 ...
    - 2.6e-8 * t_GMST^3;

THETA_GMST_deg_wo_Eq = THETA_GMST_1982_deg_0h + rad2deg(EOP_params.omega_earth_rad_sec * time_struct.sec_past_midnight_UT1);
THETA_GMST_deg_wo_Eq = mod(THETA_GMST_deg_wo_Eq,360);

% Eq 3-79:
Eq_equinox_1982_deg = ... 
                        DELTA_psi_1980_deg*cosd(epsilon_mean_1980_deg) + ...
                        0.00264*arcsec2deg*sind(Omega_deg) + ...
                        0.000063*arcsec2deg*sind(2*Omega_deg);

THETA_GAST_1982_deg = THETA_GMST_deg_wo_Eq + Eq_equinox_1982_deg;

S = rot3_deg(-THETA_GAST_1982_deg);
%% M Polar Motion Matrix

xp_rad = EOP_params.x_pole_arcsec * arcsec2rad;
yp_rad = EOP_params.y_pole_arcsec * arcsec2rad;

% Approximate Version
% M = [1 0  -xp_rad;
%      0 1    yp_rad;
%      xp_rad -yp_rad 1];

M = rot1_rad(yp_rad)*rot2_rad(xp_rad); % Full version
%% Orthogonalize All Matrices
M = Tools.orthodcm(M);
S = Tools.orthodcm(S);
N = Tools.orthodcm(N);
P = Tools.orthodcm(P);
%% Compute Transform:

Omega_Earth_Vector_Rad_Sec = [0; 0; EOP_params.omega_earth_rad_sec];

if strcmp(conversion_type,'ECEF_to_ECI') % This is ITRF to GCRF, which is the same as ECEF to ECI
    r_PEF_km = M*r_original_km;
    r_converted_km = P*N*S*M*r_original_km;
    v_converted_km_sec = P*N*S*(M*v_original_km_sec + cross(Omega_Earth_Vector_Rad_Sec,r_PEF_km));

elseif strcmp(conversion_type,'ECI_to_ECEF') % This is GCRF to ITRF, which is the same as ECI to ECEF
    r_PEF_km = S'*N'*P'*r_original_km;
    r_converted_km = M'*r_PEF_km;
    v_converted_km_sec = M'*(S'*N'*P'*v_original_km_sec - cross(Omega_Earth_Vector_Rad_Sec,r_PEF_km));
end
end
 %% Rotation Matrix Functions (Generated Quickly through Gemini formatting, but I understand their function 100%)

function R1 = rot1_deg(angle_deg)
    % Rotation about the x-axis
    theta = deg2rad(angle_deg);
    c = cos(theta);
    s = sin(theta);
    R1 = [1,  0,  0;
          0,  c,  s;
          0, -s,  c];
end

function R2 = rot2_deg(angle_deg)
    % Rotation about the y-axis
    theta = deg2rad(angle_deg);
    c = cos(theta);
    s = sin(theta);
    R2 = [ c,  0, -s;
           0,  1,  0;
           s,  0,  c];
end

function R3 = rot3_deg(angle_deg)
    % Rotation about the z-axis
    theta = deg2rad(angle_deg);
    c = cos(theta);
    s = sin(theta);
    R3 = [ c,  s,  0;
          -s,  c,  0;
           0,  0,  1];
end

function R1 = rot1_rad(angle_rad)
    % Rotation about the x-axis
    theta = angle_rad;
    c = cos(theta);
    s = sin(theta);
    R1 = [1,  0,  0;
          0,  c,  s;
          0, -s,  c];
end

function R2 = rot2_rad(angle_rad)
    % Rotation about the y-axis
    theta = angle_rad;
    c = cos(theta);
    s = sin(theta);
    R2 = [ c,  0, -s;
           0,  1,  0;
           s,  0,  c];
end

function R3 = rot3_rad(angle_rad)
    % Rotation about the z-axis
    theta = angle_rad;
    c = cos(theta);
    s = sin(theta);
    R3 = [ c,  s,  0;
          -s,  c,  0;
           0,  0,  1];
end