clear; clc;

%% SETUP -- just because we need the ENV and S stuff
case_names = {'A','B','C','D','E','F'};
scenario_config = {'Accel','Final_3D','Final_1D','HW5','24Dynamics'};
[S, ENV] = Setup.loadSettings('F','Final_1D',false,false);

%%
STM_0 = eye(6);
y0 = [S.IC_Sat_Epoch.position_ECI_meters;...
      S.IC_Sat_Epoch.velocity_ECI_meters_per_second;...
      STM_0(:)];

tic;
options = odeset('RelTol',1e-12,'AbsTol',1e-14);
tspan = 0:60:Units.SEC_IN_SOLAR_DAY;
[t, y] = ode45(@(t,X) jah_sat_1_ode( ...
        t, X, S, ENV, true), ...
        tspan, y0, options);
disp(toc)

%%
State_Ref = load("ref_data\24hourStates.mat").State * Units.KILOMETERS;
State_Isaac = y(:,1:6);

%% 
Diff_State = (State_Isaac - State_Ref) ./ State_Ref;

plot(t,Diff_State)
ylim([-1,1])



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
plot3(State_Isaac(:,1), State_Isaac(:,2), State_Isaac(:,3), 'r-', 'LineWidth', 2);
hold on;
plot3(State_Ref(:,1), State_Ref(:,2), State_Ref(:,3), 'b-', 'LineWidth', 2);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');   
axis equal;

%%
% Animate satellite positions (Isaac vs Reference) over the trajectory
figure('Name','Satellite Animation','Color','w');
% Earth
[X_e, Y_e, Z_e] = sphere(60);
R_e = PhysicsConstants.R_EARTH_KM * Units.KILOMETERS;
hSurf = surf(X_e*R_e, Y_e*R_e, Z_e*R_e, ...
    'FaceColor',[0.2 0.4 0.8], 'EdgeColor','none', 'FaceAlpha', 0.5);
hold on;
% Trajectory lines (full)
hIsaacTrail = plot3(State_Isaac(:,1), State_Isaac(:,2), State_Isaac(:,3), ...
    'r:', 'LineWidth', 1);
hRefTrail = plot3(State_Ref(:,1), State_Ref(:,2), State_Ref(:,3), ...
    'b:', 'LineWidth', 1);
% Current position markers
hIsaac = plot3(State_Isaac(1,1), State_Isaac(1,2), State_Isaac(1,3), ...
    'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
hRef = plot3(State_Ref(1,1), State_Ref(1,2), State_Ref(1,3), ...
    'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
% Lines from Earth center to satellites (optional)
hIsaacLine = plot3([0 State_Isaac(1,1)], [0 State_Isaac(1,2)], [0 State_Isaac(1,3)], ...
    'r-', 'LineWidth', 0.5);
hRefLine = plot3([0 State_Ref(1,1)], [0 State_Ref(1,2)], [0 State_Ref(1,3)], ...
    'b-', 'LineWidth', 0.5);

% Formatting
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal;
grid on;
view(3);
% Set axis limits to encompass both trajectories with margin
allPos = [State_Isaac(:,1:3); State_Ref(:,1:3)];
m = min(allPos,[],1); M = max(allPos,[],1);
pad = 0.1 * max(M - m);
xlim([m(1)-pad, M(1)+pad]); ylim([m(2)-pad, M(2)+pad]); zlim([m(3)-pad, M(3)+pad]);

% Create legend
legend([hIsaac, hRef], {'Isaac','Reference'}, 'Location','northeastoutside');

% Animation loop
nFrames = size(State_Isaac,1);
% Choose update stride to keep animation reasonable if many frames
maxFrames = 1000;
stride = max(1, ceil(nFrames / maxFrames));
for k = 1:stride:nFrames
    % Update marker positions
    set(hIsaac, 'XData', State_Isaac(k,1), 'YData', State_Isaac(k,2), 'ZData', State_Isaac(k,3));
    set(hRef,   'XData', State_Ref(k,1),   'YData', State_Ref(k,2),   'ZData', State_Ref(k,3));
    % Update radial lines
    set(hIsaacLine, 'XData', [0 State_Isaac(k,1)], 'YData', [0 State_Isaac(k,2)], 'ZData', [0 State_Isaac(k,3)]);
    set(hRefLine,   'XData', [0 State_Ref(k,1)],   'YData', [0 State_Ref(k,2)],   'ZData', [0 State_Ref(k,3)]);
    % Optionally show progress time in title
    title(sprintf('Time = %s', datestr(seconds(t(k))/86400,'HH:MM:SS')));
    drawnow;
end
