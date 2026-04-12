function P = computeLegendre(r_ecef, L_max)
    % Vallado 8-
    % Inputs:
    %   r_ecef : [x; y; z] position vector in meters
    %   L_max  : Maximum degree (e.g., 20 for EGM96 20x20)
    %
    % Output:
    %   P      : (L_max+1) x (L_max+1) matrix of Associated Legendre Functions
    
    % 1. Calculate geometry ratios to avoid trig functions
    r_mag = norm(r_ecef);
    gamma = r_ecef(3) / r_mag;         % sin(phi_gc) = z/r
    cos_phi_gc   = norm(r_ecef(1:2)) / r_mag; % cos(phi_gc) = sqrt(x^2 + y^2)/r
    
    % 2. Initialize matrix (Rows = L+1, Cols = m+1)
    P = zeros(L_max + 1, L_max + 2); ; % We let the 2nd dimension go to L_max + 2 so that we can pull out m+1 coefficients!
    
    % 3. Seed starting values (Degree 0 and 1)
    P(0+1, 0+1) = 1;                % P0,0
    P(1+1, 0+1) = gamma;            % P1,0
    P(1+1, 1+1) = cos_phi_gc;       % P1,1
    
    % 4. Recursive Computation
    for L = 2:L_max
        for m = 0:L
            % Case 1: Zonal Terms (m = 0)
            if m == 0
                P(L+1, m+1) = ((2*L - 1) * gamma * P(L-1+1, 0+1) - (L - 1) * P(L-2+1, 0+1)) / L;
            
            % Case 2: Sectorial Terms (m = L)
            elseif m == L
                P(L+1, m+1) = (2*L - 1) * cos_phi_gc * P(L-1+1, L-1+1);
            
            % Case 3: Tesseral Terms (0 < m < L)
            else
                % Using Equation 8-57: P_l,m = P_l-2,m + (2l-1)*cos(phi)*P_l-1,m-1
                P(L+1, m+1) = P(L-2+1, m+1) + (2*L - 1) * cos_phi_gc * P(L-1+1, m-1+1);
            end
        end
    end
end