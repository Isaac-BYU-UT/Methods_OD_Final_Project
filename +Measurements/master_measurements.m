function [Measurement_Computation_func, H_matrix_computation_func] = master_measurements()
%% This function computes the symbolic expressions for the range and range rate measurements, and H matrix

syms r_km [3 1] real
syms v_km_s [3 1] real
syms C_drag real
syms r_station_km [3 1] real % station position [km]
syms v_station_km [3 1] real % station velocity [km/s]

X_states = [r_km; v_km_s; C_drag];

G_Matrix = Measurements.Range_Range_Rate(r_km, v_km_s, r_station_km, v_station_km);
H_Matrix = jacobian(G_Matrix, X_states);

Measurement_Computation_func = matlabFunction(G_Matrix, ...
    'Vars', {r_km, v_km_s, r_station_km, v_station_km});

H_matrix_computation_func = matlabFunction(H_Matrix, ...
    'Vars', {r_km, v_km_s, r_station_km, v_station_km});
end