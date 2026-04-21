function [] = master_dynamics()

% Define symbolic expressions
syms r_ECI_meters [3 1] real
syms v_ECI_meters_s [3 1] real
syms C_drag real
syms r_sun_rel_earth_ECI_meters [3 1] real
syms r_moon_rel_earth_ECI_meters [3 1] real
syms sat_is_illuminated real
syms R_ECEF_from_ECI [3 3] real

% X_states = [r_ECI_meters; v_ECI_meters_s; C_drag];

% ALL ACCELERATIONS IN M/S^2 and in the ECI frame!
% ---------------------------
a_2B        = Forces.Gravity_2Body(r_ECI_meters);

a_Zonals    = Forces.Gravity_Zonal(r_ECI_meters,...
                                    R_ECEF_from_ECI,...
                                    true, true, true); % J2, J3, J4 toggles.

a_Drag      = Forces.Atmospheric_Drag(r_ECI_meters,...
                                      v_ECI_meters_s,...
                                      C_drag, ...
                                      r_sun_rel_earth_ECI_meters);

a_LuniSolar = Forces.Luni_Solar_Pertubations(r_ECI_meters, ...
                                             r_sun_rel_earth_ECI_meters,...
                                             r_moon_rel_earth_ECI_meters);

a_SRP       = Forces.Solar_Radiation_Pressure(r_ECI_meters,...
                                              r_sun_rel_earth_ECI_meters,...
                                              sat_is_illuminated);

% TODO: Add full 20x20 spherical harmonics stuff

% ----------------------------

a_component = [a_2B; a_Zonals; a_Drag; a_LuniSolar; a_SRP];
a_total_simple = a_2B + a_Zonals + a_Drag + a_LuniSolar + a_SRP; % In the future state of this, the Acceleration computation will include all the way to EGM-96 20x20, 
% Todo: Make full-fledge accelearttion model with 20x20

% F_dynamics_matrix_full = [v_ECI_meters_s; a_total_simple; 0];
% F_dynamics_matrix_simple = [v_ECI_meters_s; a_total_simple; 0];

% --- Compute Partials (The heavy lifting) ---
da_dr = jacobian(a_total_simple, r_ECI_meters);
da_dv = jacobian(a_total_simple, v_ECI_meters_s);
% da_dCd = jacobian(a_total_simple, C_drag);

% --- Generate Optimized Files ---
vars = {r_ECI_meters, v_ECI_meters_s, C_drag, ... 
        r_sun_rel_earth_ECI_meters, r_moon_rel_earth_ECI_meters, ...
        sat_is_illuminated, R_ECEF_from_ECI};

matlabFunction(da_dr, 'File', '+Forces/calc_da_dr', 'Vars', vars, 'Optimize', true);
matlabFunction(da_dv, 'File', '+Forces/calc_da_dv', 'Vars', vars, 'Optimize', true);
% matlabFunction(da_dCd, 'File', '+Forces/calc_da_dCd', 'Vars', vars, 'Optimize', true);

% Convert to matlabFunction form 
% ---------------------------------
matlabFunction(a_total_simple,'File', '+Forces/Compute_Total_Acceleration_ECI_m_s2',...
                                        'Vars', vars,"Optimize",true);

matlabFunction(a_component,'File', '+Forces/Compute_Component_Acceleration_ECI_m_s2',...
                                        'Vars', vars,"Optimize",true);

end
