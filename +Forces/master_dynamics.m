function [Acceleration_Computation_func, A_matrix_computation_func] = master_dynamics()
%% This function computes the symbolic expressions for the acceleration and A-matrix, and generates numeric functions

syms r_ECI_km [3 1] real
syms v_ECI_km_s [3 1] real
syms C_drag real
syms r_sun_rel_earth_ECI_km [3 1] real
syms r_moon_rel_earth_ECI_km [3 1] real
syms sat_is_illuminated

X_states = [r_ECI_km; v_ECI_km_s; C_drag];

a_2B = Forces.Gravity_2Body(r_ECI_km);
a_Zonals = Forces.Gravity_Zonal(r_ECI_km, true, true, true); % J2, J3, J4 toggles.
a_Drag = Forces.Atmospheric_Drag(r_ECI_km, v_ECI_km_s, C_drag, r_sun_rel_earth_ECI_km);
a_LuniSolar = Forces.Luni_Solar_Pertubations(r_ECI_km, r_sun_rel_earth_ECI_km, r_moon_rel_earth_ECI_km);
a_SRP = Forces.Solar_Radiation_Pressure(r_ECI_km, r_sun_rel_earth_ECI_km, sat_is_illuminated);

a_total = a_2B + a_Zonals + a_Drag + a_LuniSolar + a_SRP; % In the future state of this, the Acceleration computation will include all the way to EGM-96 20x20, 

F_dynamics = [v_ECI_km_s; a_total; 0];
A_matrix = jacobian(F_dynamics, X_states); % The A matrix will just use pertubations, not EGM-96 20x20.

Acceleration_Computation_func = matlabFunction(a_total, 'Vars', {r_ECI_km, v_ECI_km_s, C_drag, r_sun_rel_earth_ECI_km, r_moon_rel_earth_ECI_km,sat_is_illuminated});
A_matrix_computation_func = matlabFunction(A_matrix, 'Vars', {r_ECI_km, v_ECI_km_s, C_drag, r_sun_rel_earth_ECI_km, r_moon_rel_earth_ECI_km,sat_is_illuminated});

end
