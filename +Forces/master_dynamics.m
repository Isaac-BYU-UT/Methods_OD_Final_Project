function [Acceleration_Computation_func, A_matrix_computation_func] = master_dynamics()
%% This function computes the symbolic expressions for the acceleration and A-matrix, and generates numeric functions

syms r_km [3 1] real
syms v_km_s [3 1] real
syms C_drag real

X_states = [r_km; v_km_s; C_drag];

a_2B = Forces.Gravity_2Body(r_km);

a_Zonals = Forces.Gravity_Zonal(r_km, true, false, false); % J2, J3, J4 toggles.

a_Drag = Forces.Atmospheric_Drag(r_km, v_km_s, C_drag);

a_total = a_2B + a_Zonals + a_Drag;

F_dynamics = [v_km_s; a_total; 0];

A_matrix = jacobian(F_dynamics, X_states);

% In the future state of this, the Acceleration computation will include all the way to EGM-96 20x20, but the A matrix will just use the simpler Jacobians

Acceleration_Computation_func = matlabFunction(a_total, 'Vars', {r_km, v_km_s, C_drag});
A_matrix_computation_func = matlabFunction(A_matrix, 'Vars', {r_km, v_km_s, C_drag});

end
