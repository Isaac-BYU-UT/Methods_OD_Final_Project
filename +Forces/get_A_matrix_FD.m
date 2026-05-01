function A = get_A_matrix_FD(t, X, S, ENV, h)
    % GET_A_MATRIX_FD Computes the A matrix (state transition matrix Jacobian) 
    %                  using finite differences instead of symbolic computation.
    %
    % This avoids the 20+ minute symbolic compilation time and is more robust
    % for complex dynamics models.
    %
    % Inputs:
    %   t      - Current time (seconds from epoch)
    %   X      - State vector [r_ECI (m), v_ECI (m/s), STM (36 elems), ...]
    %   S      - Settings structure
    %   ENV    - Environment structure
    %   h      - Finite difference step size (optional, default 1e-6)
    %
    % Output:
    %   A      - 6x6 state transition matrix Jacobian
    %            [  0  I  ]
    %            [ da/dr da/dv ]
    %
    % Usage:
    %   In jah_sat_1_ode.m:
    %   A_matrix = Forces.get_A_matrix_FD(t, X, S, ENV);
    %   STM_dot = A_matrix * STM;

    if nargin < 5
        h = 1e-6;  % Default step size
    end
    
    % Extract current state
    r_ECI = X(1:3);
    v_ECI = X(4:6);
    
    % Initialize A matrix (6x6)
    A = zeros(6, 6);
    
    % ---- Upper right block: dr/dv = Identity ----
    A(1:3, 4:6) = eye(3);
    
    % ---- Compute da/dr using central differences ----
    for i = 1:3
        X_plus = X;
        X_minus = X;
        X_plus(i) = r_ECI(i) + h;
        X_minus(i) = r_ECI(i) - h;
        
        a_plus = Forces.compute_total_acceleration(t, X_plus, S, ENV);
        a_minus = Forces.compute_total_acceleration(t, X_minus, S, ENV);
        
        A(4:6, i) = (a_plus - a_minus) / (2 * h);
    end
    
    % ---- Compute da/dv using central differences ----
    for i = 4:6
        X_plus = X;
        X_minus = X;
        X_plus(i) = v_ECI(i-3) + h;
        X_minus(i) = v_ECI(i-3) - h;
        
        a_plus = Forces.compute_total_acceleration(t, X_plus, S, ENV);
        a_minus = Forces.compute_total_acceleration(t, X_minus, S, ENV);
        
        A(4:6, i-3) = (a_plus - a_minus) / (2 * h);
    end
    
end
