function plot_RIC_ellipses_all_cases(mat_path)

clc; close all;

%% Load data
data = load(mat_path);

cases = {'A','B','C','D','E','F','G'};

for k = 1:length(cases)
    c = cases{k};
    pos.(c) = data.(['sorensen_pos_case' c]);
    vel.(c) = data.(['sorensen_vel_case' c]);
    P.(c)   = data.(['sorensen_poscov_case' c]);
end

%% Truth (Case F)
r_ref = pos.F;
v_ref = vel.F;

if any(isnan(r_ref)) || any(isnan(v_ref))
    error('Case F missing');
end

%% Build RIC frame
R_hat = r_ref / norm(r_ref);
h_vec = cross(r_ref, v_ref);
C_hat = h_vec / norm(h_vec);
I_hat = cross(C_hat, R_hat);

T = [R_hat'; I_hat'; C_hat'];  % ECI → RIC

planes = {
    [1 2], 'R-I';
    [1 3], 'R-C';
    [2 3], 'I-C'
};

colors = lines(length(cases));

%% Loop over planes → ONE FIGURE EACH
for p = 1:3
    
    idx = planes{p,1};
    name = planes{p,2};
    
    figure; hold on; grid on;
    title([name ' Plane (RIC)']);
    xlabel(name(1)); ylabel(name(3));
    
    for k = 1:length(cases)
        c = cases{k};
        
        r = pos.(c);
        P_case = P.(c);
        
        if any(isnan(r)) || any(isnan(P_case(:)))
            continue;
        end
        
        %% Position difference (relative to F)
        dr = r - r_ref;
        dr_RIC = T * dr;
        
        %% Covariance transform
        P_rr = P_case(1:3,1:3);
        P_RIC = T * P_rr * T';
        
        %% Extract 2D slice
        mu = dr_RIC(idx);
        P2 = P_RIC(idx, idx);
        
        %% Plot ellipse (3σ)
        plot_cov_ellipse(mu, P2, 3, colors(k,:));
        
        %% Plot center
        plot(mu(1), mu(2), 'o', 'Color', colors(k,:), ...
            'DisplayName', ['Case ' c]);
    end
    
    %% Plot truth at origin
    plot(0,0,'kp','MarkerSize',10,'MarkerFaceColor','k', ...
        'DisplayName','Case F (Truth)');
    
    
    legend show;

end

end


function plot_cov_ellipse(mu, P, nsig, color)

theta = linspace(0, 2*pi, 100);

[V, D] = eig(P);
A = nsig * V * sqrt(D);

ellipse = A * [cos(theta); sin(theta)];

x = ellipse(1,:) + mu(1);
y = ellipse(2,:) + mu(2);

plot(x, y, 'Color', color, 'LineWidth', 1.5);

end