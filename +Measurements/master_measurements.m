function [] = master_measurements()
%% This function computes the symbolic expressions for the range and range rate measurements, and H matrix

syms r_ECI_meters [3 1] real
syms v_ECI_meters_s [3 1] real
% syms C_drag real
syms r_station_ECI_meters [3 1] real % station position [m]
syms v_station_ECI_meters [3 1] real % station velocity [m/s]

X_states = [r_ECI_meters; v_ECI_meters_s]; % Removed drag from state

G_Matrix = Measurements.Compute_Range_Range_Rate(r_ECI_meters, v_ECI_meters_s, r_station_ECI_meters, v_station_ECI_meters);
H_Matrix = jacobian(G_Matrix, X_states);

matlabFunction(H_Matrix, ...
    'File', '+Measurements/Compute_H_matrix',...
    'Vars', {r_ECI_meters, v_ECI_meters_s, r_station_ECI_meters, v_station_ECI_meters},...
    'Optimize',true);
end