clear; clc;

%% SETUP
case_names = {'A','B','C','D','E','F'};
scenario_config = {'Accel','Final_3D','Final_1D','HW5'};
[S, ENV] = Setup.loadSettings('F','Final_1D',false,false);

%%

initial_state = EKF.initialize_ekf(S, ENV);
tic;
options = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t1, y1] = ode113(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, false), ...
        [0, 5000], initial_state.X_input, options);
t1 = toc;
disp(t1);

options = odeset('RelTol',1e-9,'AbsTol',1e-11);
[t, y] = ode45(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, false), ...
        [0, 5000], initial_state.X_input, options);
disp(toc)

y1(end,1:3) - y(end,1:3)

%% Plot

% Plot the orbital trajectory with the earth
figure;
hold on;
% Earth
[X_e, Y_e, Z_e] = sphere(50);
R_e = PhysicsConstants.R_EARTH_KM * Units.KILOMETERS;
surf(X_e*R_e, Y_e*R_e, Z_e*R_e, ...
    'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% Satellite trajectory
plot3(y(:,1), y(:,2), y(:,3), 'r-', 'LineWidth', 2);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');   
axis equal;

