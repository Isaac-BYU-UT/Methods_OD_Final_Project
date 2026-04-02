function plot_position(r_ECI, v_ECI, P_ECI, full_orbit_ECI, orbit_on)
    % r_ECI, v_ECI: Current state (column vectors)
    % P_ECI: Covariance
    % full_orbit_ECI: [Nx3] matrix of the orbit path in ECI
    
    r_ECI = r_ECI(:); v_ECI = v_ECI(:);
    R_earth = Constants.R_EARTH_KM;

    % 1. Define RSW Basis and Rotation
    r_unit = r_ECI / norm(r_ECI);
    w_unit = cross(r_ECI, v_ECI) / norm(cross(r_ECI, v_ECI));
    s_unit = cross(w_unit, r_unit);
    T = [r_unit'; s_unit'; w_unit']; % ECI -> RSW

    % 2. Transform Earth and Orbit into RSW
    % Earth Center in RSW
    earth_rsw = T * (-r_ECI); 
    
    % Orbit Path in RSW (Relative to current satellite position)
    orbit_rsw = (T * (full_orbit_ECI(:,1:3) - r_ECI')')';

    % 3. Covariance
    P_RSW = T * P_ECI(1:3, 1:3) * T';

    % 4. Plotting
    labels = {'Radial (R)', 'Intrack (S)', 'Crosstrack (W)'};
    pairs = [1,2; 1,3; 3,2]; 
    titles = {'Radial vs Intrack', 'Radial vs Crosstrack', 'Crosstrack vs Intrack'};
    
    figure('Color', 'w');
    
    for i = 1:3
        subplot(1, 3, i); hold on;
        idx = pairs(i,:);
        
        % Plot Earth Disk (Approximate)
        t_circ = linspace(0, 2*pi, 100);
        fill(earth_rsw(idx(1)) + R_earth*cos(t_circ), ...
             earth_rsw(idx(2)) + R_earth*sin(t_circ), ...
             [0.8 0.9 1], 'EdgeColor', [0.2 0.4 0.8], 'DisplayName', 'Earth');

        % Plot Orbit Path
        if orbit_on
            plot(orbit_rsw(:,idx(1)), orbit_rsw(:,idx(2)), 'k--', 'LineWidth', 1, 'DisplayName', 'Orbit');
        end
        
        % Plot 3-Sigma Covariance (The "Useful" part)
        error_ellipse(P_RSW(idx, idx), [0,0], 3);
        
        xlabel(labels{idx(1)}); ylabel(labels{idx(2)});
        title(titles{i});
        grid on; axis equal;
        
        % Zoom in on the satellite/covariance
        % Adjust this multiplier to see more or less of the orbit
        limit = 50; % km
        xlim([-limit, limit]); ylim([-limit, limit]);
    end
end

function error_ellipse(C, mu, n_sigma)
    [V, D] = eig(C);
    t = linspace(0, 2*pi, 100);
    
    % Ensure eigenvalues are positive (numerical precision check)
    eigs = max(diag(D), 0);
    
    % Create circle and scale by n_sigma * sigma
    xy = [sqrt(eigs(1))*cos(t); sqrt(eigs(2))*sin(t)] * n_sigma;
    
    % Rotate to ECI/RSW orientation
    ellipse = V * xy + mu(:);
    
    plot(ellipse(1,:), ellipse(2,:), 'LineWidth', 2, 'Color', [0 .45 .74]);
    hold on;
    plot(mu(1), mu(2), 'r+', 'MarkerSize', 10); 
end