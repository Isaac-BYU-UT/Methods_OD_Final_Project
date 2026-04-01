function plot_orbit_errors(r_ECI, v_ECI, P_ECI)
    % r_ECI: [3x1] Position (km)
    % v_ECI: [3x1] Velocity (km/s)
    % P_ECI: [6x6] or [3x3] Covariance matrix (km^2)
    
    % 1. Define RSW Basis Vectors
    r_unit = r_ECI / norm(r_ECI);
    h_vec  = cross(r_ECI, v_ECI);
    w_unit = h_vec / norm(h_vec);
    s_unit = cross(w_unit, r_unit);
    
    % 2. Construct Rotation Matrix (ECI to RSW)
    % Note: P is 6x6, but we usually only plot position errors (top-left 3x3)
    T = [r_unit'; s_unit'; w_unit'];
    P_pos_ECI = P_ECI(1:3, 1:3);
    P_RSW = T * P_pos_ECI * T';
    
    % 3. Setup Plotting
    labels = {'Radial (km)', 'Intrack (km)', 'Crosstrack (km)'};
    pairs = [1,2; 1,3; 3,2]; % R-S, R-W, W-S
    titles = {'Radial vs Intrack', 'Radial vs Crosstrack', 'Crosstrack vs Intrack'};
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 400]);
    
    for i = 1:3
        subplot(1, 3, i);
        idx = pairs(i,:);
        
        % Extract 2x2 sub-covariance for the plane
        sub_P = P_RSW(idx, idx);
        
        % Generate Ellipse Points (3-sigma)
        error_ellipse(sub_P, [0,0], 3);
        
        xlabel(labels{idx(1)});
        ylabel(labels{idx(2)});
        title(titles{i});
        grid on; axis equal;
    end
end

function error_ellipse(C, mu, n_sigma)
    % Helper to plot a 2D ellipse
    [V, D] = eig(C);
    t = linspace(0, 2*pi, 100);
    % Scaled circle
    xy = [cos(t); sin(t)] * n_sigma;
    % Transform to ellipse: rotate by eigenvectors, scale by sqrt(eigenvalues)
    ellipse = V * sqrt(D) * xy + mu';
    
    plot(ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', [0 .45 .74]);
    hold on;
    plot(mu(1), mu(2), 'r+', 'MarkerSize', 10); % Mean
end