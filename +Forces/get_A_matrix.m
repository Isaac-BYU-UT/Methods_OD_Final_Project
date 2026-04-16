function A = get_A_matrix(r, v, Cd, r_sun, r_moon, illuminated, R_ECEF)
    % Initialize Sparse Matrix (7x7)
    A = zeros(7, 7);
    
    % Upper Right: dr/dv = Identity
    A(1:3, 4:6) = eye(3);
    
    % Middle Row: Partial derivatives of acceleration
    % These call the optimized .m files generated in Step 1
    A(4:6, 1:3) = Forces.calc_da_dr(r, v, Cd, r_sun, r_moon, illuminated, R_ECEF);
    A(4:6, 4:6) = Forces.calc_da_dv(r, v, Cd, r_sun, r_moon, illuminated, R_ECEF);
    A(4:6, 7)   = Forces.calc_da_dCd(r, v, Cd, r_sun, r_moon, illuminated, R_ECEF);
    
    % Bottom Row: d(Cd)/dt = 0 (Already zeroed by initialization)
end