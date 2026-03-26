function dXdt = jah_sat_1_ode(~, X, accel_func,A_func)

    r = X(1:3);
    v = X(4:6);
    C_drag = X(7);
    STM = reshape(X(8:56),7,7); % 7 states, therefor this will be a 7x7 STM

    A = A_func(r,v,C_drag);

    STM_dot = A*STM; % Propogate STM Matrix
    
    a_total = accel_func(r,v,C_drag);

    C_drag_dot = 0; % The time derivative of the drag coefficient should be 0.

    dXdt = [v; a_total; C_drag_dot; STM_dot(:)];
end