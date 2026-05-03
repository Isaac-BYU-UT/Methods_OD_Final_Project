function plot3DPositions_multi(matFileNames)
    % plot3DPositions - Plots ECI positions from multiple .mat files
    % Input: matFileNames - A cell array of strings, e.g., {'file1.mat', 'file2.mat'}
    
    if ischar(matFileNames) || isstring(matFileNames)
        matFileNames = {matFileNames}; % Ensure it's a cell array for iteration
    end

    % Constants
    earth_radius = 6378137; % Earth radius in meters
    
    % Setup Figure
    figure('Color', 'w');
    hold on; grid on; axis equal;
    
    % Get default MATLAB color order to cycle through colors automatically
    colors = get(gca, 'ColorOrder');
    numColors = size(colors, 1);
    
    legendEntries = []; % To store plot handles for the legend

    % Iterate through each file
    for f = 1:length(matFileNames)
        currentFile = matFileNames{f};
        data = load(currentFile);
        vars = fieldnames(data);
        
        x_coords = [];
        y_coords = [];
        z_coords = [];
        labels = {};
        
        % Select color (cycle if more files than colors)
        colorIdx = mod(f-1, numColors) + 1;
        currentColor = colors(colorIdx, :);

        % Extract variables
        for i = 1:length(vars)
            if contains(vars{i}, 'sorensen_pos_case') && ~contains(vars{i}, 'poscov')
                pos_vector = data.(vars{i});
                x_coords(end+1) = pos_vector(1);
                y_coords(end+1) = pos_vector(2);
                z_coords(end+1) = pos_vector(3);
                
                % Clean up label for individual points
                labels{end+1} = strrep(vars{i}, 'sorensen_pos_', '');
            end
        end
        
        % Plot the set of positions for this file
        h = scatter3(x_coords, y_coords, z_coords, 100, currentColor, 'filled', ...
                     'MarkerEdgeColor', 'k', 'DisplayName', currentFile);
        legendEntries(end+1) = h; % Store handle for legend
        
        % Add Labels to points with a slight offset
        text(x_coords, y_coords, z_coords, labels, 'VerticalAlignment', 'bottom', ...
             'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', currentColor * 0.7);
    end

    % % --- Optional: Plot Earth Wireframe ---
    % [sx, sy, sz] = sphere(50);
    % mesh(sx*earth_radius, sy*earth_radius, sz*earth_radius, ...
    %      'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'HandleVisibility', 'off');

    % Formatting
    xlabel('X (ECI) [m]');
    ylabel('Y (ECI) [m]');
    zlabel('Z (ECI) [m]');
    title('Multi-File Sorensen Satellite Positions');
    view(3);
    
    % Add legend using the file names
    legend(legendEntries, 'Interpreter', 'none', 'Location', 'bestoutside');
end