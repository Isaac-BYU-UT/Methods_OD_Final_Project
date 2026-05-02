function X_state_dot = jah_sat_1_ode_OptionB(t, X, S, ENV, debug_on)

    if nargin < 4
        debug_on = false;
    end
    
    r_ECI_meters = X(1:3);
    v_ECI_meters_s = X(4:6);
    C_drag = SatelliteProperties.C_Drag;
    STM = reshape(X(7:42),6,6); % Now only 6 states! % 7 states, therefor this will be a 7x7 STM


    % Compute and Propogate STM
    % -------------------------
    % Using finite difference A-matrix to avoid 20+ min symbolic compilation
    A_matrix = Forces.get_A_matrix_FD(t, X, S, ENV);
    
    STM_dot = A_matrix*STM; 

    if t == 0
        disp("A_matrix: "); disp((A_matrix - S.ref_data.A_t0_ref(1:6,1:6))./S.ref_data.A_t0_ref(1:6,1:6));
    end

    a_total_ECI_m_s2 = Forces.compute_total_acceleration(t, X, S, ENV);

    % C_drag_dot = 0;

    X_state_dot =  [v_ECI_meters_s;...
                    a_total_ECI_m_s2;...
                    % C_drag_dot;...
                    STM_dot(:)];



end
