function plot3DPositions(matFileName)
    % Load the .mat file
    data = load(matFileName);
    vars = fieldnames(data);
    
    % Constants
    earth_radius = 6378137; % Earth radius in meters
    
    % Initialize arrays for plotting
    x_coords = [];
    y_coords = [];
    z_coords = [];
    labels = {};

    % Extract variables matching the "sorensen_pos_case" pattern
    for i = 1:length(vars)
        if contains(vars{i}, 'sorensen_pos_case') && ~contains(vars{i}, 'poscov')
            pos_vector = data.(vars{i});
            x_coords(end+1) = pos_vector(1);
            y_coords(end+1) = pos_vector(2);
            z_coords(end+1) = pos_vector(3);
            % Clean up label (e.g., 'caseA')
            labels{end+1} = strrep(vars{i}, 'sorensen_pos_', '');
        end
    end

    % --- Plotting ---
    figure('Color', 'w');
    hold on; grid on; axis equal;
    
    % % 1. Plot Earth (Wireframe/Surface)
    % [sx, sy, sz] = sphere(50);
    % surf(sx*earth_radius, sy*earth_radius, sz*earth_radius, ...
    %     'FaceColor', 'blue', 'EdgeColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1);

    % 2. Plot Positions
    scatter3(x_coords, y_coords, z_coords, 100, 'red', 'filled', 'MarkerEdgeColor', 'k');
    
    % 3. Add Labels to points
    text(x_coords, y_coords, z_coords, labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    % Formatting
    xlabel('X (ECI) [m]');
    ylabel('Y (ECI) [m]');
    zlabel('Z (ECI) [m]');
    title('Sorensen Satellite Positions in ECI Frame');
    view(3); % Set 3D view
    
    % Add a legend for clarity
    legend('Earth', 'Case Positions', 'Location', 'bestoutside');
end