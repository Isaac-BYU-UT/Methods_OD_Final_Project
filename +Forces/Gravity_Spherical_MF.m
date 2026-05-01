function accel_ECI_spherical_harmonics_m_s2 = Gravity_Spherical_MF(r_ECI_meters, ...
                                                              R_ECEF_from_ECI, ...
                                                              C_coeff_matrix, ...
                                                              S_coeff_matrix, ...
                                                              P_coeff_matrix)
% Gravity_Spherical_MF
% A pure, MatlabFunction-friendly spherical harmonics gravity model.
% This version is self-contained and avoids external state such as IC/ENV.

%% Constants
L_max = 20;                         % 20x20 EGM96
mu = PhysicsConstants.MU_EARTH_KM3_S2;
R_earth_meters = PhysicsConstants.R_EARTH_KM * Units.KILOMETERS;

%% Position in ECEF
r_ECEF_meters = R_ECEF_from_ECI * r_ECI_meters;
r = norm(r_ECEF_meters);

%% Longitude / latitude preparation
x_y_norm_meters = sqrt(r_ECEF_meters(1)^2 + r_ECEF_meters(2)^2);
cos_lambda = r_ECEF_meters(1) / x_y_norm_meters;
sin_lambda = r_ECEF_meters(2) / x_y_norm_meters;
tan_phi_gc = r_ECEF_meters(3) / x_y_norm_meters;


if isa(cos_lambda, 'sym')
    Cos_m_lambda = sym(zeros(1, L_max + 1));
    Sin_m_lambda = sym(zeros(1, L_max + 1));
    m_Tan_phi_gc = sym(zeros(1, L_max + 1));
else
    Cos_m_lambda = zeros(1, L_max + 1);
    Sin_m_lambda = zeros(1, L_max + 1);
    m_Tan_phi_gc = zeros(1, L_max + 1);
end

Cos_m_lambda(1) = 1;
Sin_m_lambda(1) = 0;
m_Tan_phi_gc(1) = 0;

if L_max >= 1
    Cos_m_lambda(2) = cos_lambda;
    Sin_m_lambda(2) = sin_lambda;
    m_Tan_phi_gc(2) = tan_phi_gc;
end

for m = 3:(L_max + 1)
    Cos_m_lambda(m) = 2 * cos_lambda * Cos_m_lambda(m-1) - Cos_m_lambda(m-2);
    Sin_m_lambda(m) = 2 * cos_lambda * Sin_m_lambda(m-1) - Sin_m_lambda(m-2);
    m_Tan_phi_gc(m) = (m-1) * tan_phi_gc + tan_phi_gc;
end

%% Partials initialization

dU_dr = 0;
dU_dphi = 0;
dU_dlambda = 0;

rI = r_ECEF_meters(1);
rJ = r_ECEF_meters(2);
rK = r_ECEF_meters(3);
rho_sq = rI^2 + rJ^2;

for L = 2:L_max
    r_ratio_pow = (R_earth_meters / r)^L;
    L_idx = L + 1;

    for m = 0:L
        m_idx = m + 1;

        P_lm = P_coeff_matrix(L_idx, m_idx);
        C_lm = C_coeff_matrix(L_idx, m_idx);
        S_lm = S_coeff_matrix(L_idx, m_idx);

        cos_mlam = Cos_m_lambda(m_idx);
        sin_mlam = Sin_m_lambda(m_idx);
        trig_sum = C_lm * cos_mlam + S_lm * sin_mlam;

        dU_dr = dU_dr + (r_ratio_pow / r) * (L + 1) * P_lm * trig_sum;

        if m < L
            P_l_mplus1 = P_coeff_matrix(L_idx, m_idx + 1);
        else
            P_l_mplus1 = 0;
        end
        m_tan_phi = m_Tan_phi_gc(m_idx);
        dU_dphi = dU_dphi + r_ratio_pow * (P_l_mplus1 - m_tan_phi * P_lm) * trig_sum;

        trig_diff = S_lm * cos_mlam - C_lm * sin_mlam;
        dU_dlambda = dU_dlambda + r_ratio_pow * m * P_lm * trig_diff;
    end
end

%% Apply constants from Vallado Eq. 8-25

dU_dr      = -(mu / r^2) * dU_dr;
dU_dphi    = (mu / r) * dU_dphi;
dU_dlambda = (mu / r) * dU_dlambda;

%% ECEF acceleration from spherical harmonic partials

a_I = ((1/r * dU_dr) - (rK / (r^2 * sqrt(rho_sq)) * dU_dphi)) * rI - (dU_dlambda / rho_sq) * rJ - (mu * rI / r^3);
a_J = ((1/r * dU_dr) - (rK / (r^2 * sqrt(rho_sq)) * dU_dphi)) * rJ + (dU_dlambda / rho_sq) * rI - (mu * rJ / r^3);
a_K = ((1/r * dU_dr) * rK) + (sqrt(rho_sq) / r^2 * dU_dphi) - (mu * rK / r^3);

a_ECEF = [a_I; a_J; a_K];

%% Rotate back to ECI using the transpose of the input matrix
accel_ECI_spherical_harmonics_m_s2 = R_ECEF_from_ECI' * a_ECEF;
end