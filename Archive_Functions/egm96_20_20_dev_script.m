clear; clc;

IC = Initial_Conditions.initial_conditions_Phase1();
ENV = Environment.setup_environment_full_1Day_Subset(IC);

[IC.r_sat_ECEF_km, IC.v_sat_ECEF_km_s] = Tools.ECEF_ECI_Converter(IC.r_sat_ECI_km,IC.v_sat_ECI_km_s,IC.UTC_epoch,"ECI_to_ECEF",ENV.EOP_t0);

r_ECEF_km = IC.r_sat_ECEF_km;
v_ECEF_km_s = IC.v_sat_ECEF_km_s;

U_total = 0;
%% Constants
L_max = 20; % For 20x20 EGM 96.
mu = Constants.MU_EARTH_KM3_S2;
R_earth = Constants.R_EARTH_KM;
r_mag_ECEF_km = norm(r_ECEF_km);
r = r_mag_ECEF_km;
P_coeff_matrix = computeLegendre(r_ECEF_km, L_max);


%%
    Cos_m_lambda = zeros(1, L_max + 1);
    Sin_m_lambda = zeros(1, L_max + 1);
    m_Tan_phi_gc = zeros(1, L_max + 1);

    cos_lambda = r_ECEF_km(1)/norm(r_ECEF_km(1:2));
    sin_lambda = r_ECEF_km(2)/norm(r_ECEF_km(1:2));
    tan_phi_gc = r_ECEF_km(3)/norm(r_ECEF_km(1:2));
    
    % Base cases
    Cos_m_lambda(1) = 1; % m = 0: cos(0)
    Sin_m_lambda(1) = 0; % m = 0: sin(0)
    m_Tan_phi_gc(1) = 0; % m = 0, 0*tan(phi_gc)
    
    if L_max >= 1
        Cos_m_lambda(2) = cos_lambda; % m = 1 cos(lambda), longitude angle x/(hypotenuse)
        Sin_m_lambda(2) = sin_lambda; % m = 1 sin(lambda), longitude angle y/(hypotenuse)
        m_Tan_phi_gc(2) = tan_phi_gc; % m = 1 tan(phi_gc)
    end
    
    % Recursion for m >= 2 (Eqn. 8-58)
    for m = 3 : (L_max + 1)
        % m_idx-1 corresponds to the 'm' in the formula
        Cos_m_lambda(m) = 2 * cos_lambda * Cos_m_lambda(m-1) - Cos_m_lambda(m-2);
        Sin_m_lambda(m) = 2 * cos_lambda * Sin_m_lambda(m-1) - Sin_m_lambda(m-2);
        m_Tan_phi_gc(m) = (m-1) * tan_phi_gc + tan_phi_gc;
    end

 %%

sum_grav_potential = 0;
    
    for L = 2:L_max
        r_ratio_pow = (R_earth / r)^L;
        inner_m_sum = 0;
        
        % Loop through order m (0 to l)
        for m = 0:L
            % indices are offset by 1 because MATLAB is 1-indexed
            L_idx = L + 1;
            m_idx = m + 1;
            
            % Retrieve P_lm, C_lm, and S_lm from your matrices
            P_lm = P_coeff_matrix(L_idx, m_idx);
            C_lm = C_coeff_matrix(L_idx, m_idx);
            S_lm = S_coeff_matrix(L_idx, m_idx);
            
            % Calculate the term for this m
            inner_m_sum = inner_m_sum + P_lm * (C_lm * Cos_m_lambda(m_idx) + S_lm * Sin_m_lambda(m_idx));
        end
        
        % Add the contribution of degree l to the total
        sum_grav_potential = sum_grav_potential + r_ratio_pow * inner_m_sum;
    end
    
    % Final Potential Equation
    U = (mu / r) * (1 + sum_grav_potential);

 %% Partials Initialization (pg. 550 vallado). The rest here compiled with help from Gemini.
dU_dr = 0;
dU_dphi = 0;
dU_dlambda = 0;

% Position components for Eq (8-27)
rI = r_ECEF_km(1);
rJ = r_ECEF_km(2);
rK = r_ECEF_km(3);
rho_sq = rI^2 + rJ^2; % Square of the distance in the Equatorial plane

for L = 2:L_max
    r_ratio_pow = (R_earth / r)^L;
    
    for m = 0:L
        L_idx = L + 1;
        m_idx = m + 1;
        
        % Data retrieval
        P_lm = P_coeff_matrix(L_idx, m_idx);
        C_lm = C_coeff_matrix(L_idx, m_idx);
        S_lm = S_coeff_matrix(L_idx, m_idx);
        
        % Common Trig Terms
        cos_mlam = Cos_m_lambda(m_idx);
        sin_mlam = Sin_m_lambda(m_idx);
        trig_sum = C_lm * cos_mlam + S_lm * sin_mlam;
        
        %% 1. Radial Partial (Eq 8-25)
        dU_dr = dU_dr + (r_ratio_pow / r) * (L + 1) * P_lm * trig_sum; % Will be zeroed out if m > L
        
        %% 2. Latitudinal Partial (Eq 8-25)
        % Note: Requires P(l, m+1). Ensure your computeLegendre handles m+1
        P_l_mplus1 = P_coeff_matrix(L_idx, m_idx + 1); 
        m_tan_phi = m_Tan_phi_gc(m_idx); 
        
        dU_dphi = dU_dphi + r_ratio_pow * (P_l_mplus1 - m_tan_phi * P_lm) * trig_sum; % Will be zeroed out if m > L
        
        %% 3. Longitudinal Partial (Eq 8-25)
        trig_diff = S_lm * cos_mlam - C_lm * sin_mlam;
        dU_dlambda = dU_dlambda + r_ratio_pow * m * P_lm * trig_diff; % Will be zeroed out if m > L
    end
end

% Apply the constants leading the summations in Eq (8-25)
dU_dr      = -(mu / r^2) * dU_dr;
dU_dphi    = (mu / r) * dU_dphi;
dU_dlambda = (mu / r) * dU_dlambda;

%% 4. Combine into Cartesian Acceleration (Eq 8-27)
% This is the acceleration in the ECEF frame
a_I = ( (1/r * dU_dr) - (rK / (r^2 * sqrt(rho_sq)) * dU_dphi) ) * rI - (dU_dlambda / rho_sq) * rJ - (mu * rI / r^3);
a_J = ( (1/r * dU_dr) - (rK / (r^2 * sqrt(rho_sq)) * dU_dphi) ) * rJ + (dU_dlambda / rho_sq) * rI - (mu * rJ / r^3);
a_K = ( (1/r * dU_dr) * rK ) + (sqrt(rho_sq) / r^2 * dU_dphi) - (mu * rK / r^3);

a_ECEF = [a_I; a_J; a_K];

%% 5. Rotate back to ECI
% Use your existing converter logic to rotate the vector
% We only need the rotation matrix R_ecef2eci
R_ecef2eci = Tools.R_ECEF_from_ECI_Matrix(IC.UTC_epoch, ENV.EOP_t0);
a_ECI = R_ecef2eci * a_ECEF;


[gx,gy,gz] = gravitysphericalharmonic(...
    transpose(r_ECEF_km * 1000),'EGM96',20,'Warning');
