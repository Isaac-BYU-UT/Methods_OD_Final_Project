function accel_ECI_spherical_harmonics_m_s2 = Gravity_Spherical(r_ECI_meters,...
                                                                R_ECEF_from_ECI,...
                                                                C_coeff_matrix,...
                                                                S_coeff_matrix,...
                                                                P_coeff_matrix)

%% Initial Position in ECEF
r_ECEF_meters = R_ECEF_from_ECI * r_ECI_meters;

%% Constants
L_max = 20; % For 20x20 EGM 96.
mu = PhysicsConstants.MU_EARTH_KM3_S2 * (Units.KILOMETERS ^ 3); % Get's into m^3/s^2
R_earth_meters = PhysicsConstants.R_EARTH_KM * Units.KILOMETERS;
r_mag_ECEF_meters = norm(r_ECEF_meters);
r = r_mag_ECEF_meters;

%%
    Cos_m_lambda = zeros(1, L_max + 1);
    Sin_m_lambda = zeros(1, L_max + 1);
    m_Tan_phi_gc = zeros(1, L_max + 1);
    
    x_y_norm_meters = sqrt(r_ECEF_meters(1)^2 + r_ECEF_meters(2)^2);
    cos_lambda = r_ECEF_meters(1)/x_y_norm_meters;
    sin_lambda = r_ECEF_meters(2)/x_y_norm_meters;
    tan_phi_gc = r_ECEF_meters(3)/x_y_norm_meters;
    
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


 %% Partials Initialization (pg. 550 vallado). The rest here compiled with help from Gemini.
dU_dr = 0;
dU_dphi = 0;
dU_dlambda = 0;

% Position components for Eq (8-27)
rI = r_ECEF_meters(1);
rJ = r_ECEF_meters(2);
rK = r_ECEF_meters(3);
rho_sq = rI^2 + rJ^2; % Square of the distance in the Equatorial plane

for L = 2:L_max
    r_ratio_pow = (R_earth_meters / r)^L;
    
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
% R_ecef2eci = Tools.R_ECEF_from_ECI_Matrix(IC.UTC_epoch, ENV.EOP_t0);
accel_ECI_spherical_harmonics_m_s2 = transpose(R_ECEF_from_ECI) * a_ECEF;


end

function P = computeLegendre(r_ecef, L_max)
    % Vallado 8-
    % Inputs:
    %   r_ecef : [x; y; z] position vector in meters
    %   L_max  : Maximum degree (e.g., 20 for EGM96 20x20)
    %
    % Output:
    %   P      : (L_max+1) x (L_max+1) matrix of Associated Legendre Functions

    % 1. Calculate geometry ratios to avoid trig functions
    r_mag = norm(r_ecef);
    gamma = r_ecef(3) / r_mag;         % sin(phi_gc) = z/r
    cos_phi_gc   = norm(r_ecef(1:2)) / r_mag; % cos(phi_gc) = sqrt(x^2 + y^2)/r

    % 2. Initialize matrix (Rows = L+1, Cols = m+1)
    P = zeros(L_max + 1, L_max + 2); % We let the 2nd dimension go to L_max + 2 so that we can pull out m+1 coefficients!

    % 3. Seed starting values (Degree 0 and 1)
    P(0+1, 0+1) = 1;                % P0,0
    P(1+1, 0+1) = gamma;            % P1,0
    P(1+1, 1+1) = cos_phi_gc;       % P1,1

    % 4. Recursive Computation
    for L = 2:L_max
        for m = 0:L
            % Case 1: Zonal Terms (m = 0)
            if m == 0
                P(L+1, m+1) = ((2*L - 1) * gamma * P(L-1+1, 0+1) - (L - 1) * P(L-2+1, 0+1)) / L;

            % Case 2: Sectorial Terms (m = L)
            elseif m == L
                P(L+1, m+1) = (2*L - 1) * cos_phi_gc * P(L-1+1, L-1+1);

            % Case 3: Tesseral Terms (0 < m < L)
            else
                % Using Equation 8-57: P_l,m = P_l-2,m + (2l-1)*cos(phi)*P_l-1,m-1
                P(L+1, m+1) = P(L-2+1, m+1) + (2*L - 1) * cos_phi_gc * P(L-1+1, m-1+1);
            end
        end
    end
end