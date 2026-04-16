function accel_ECI_meters_s2 = Gravity_Zonal(r_ECI_meters, R_ECEF_from_ECI, J2_on, J3_on, J4_on)
    % GRAVITY_ZONAL Calculates acceleration due to J2, J3, and J4 zonals.
    % Refactored to match Vallado/Escobal algebraic forms.
    
    % 1. Constants
    mu = PhysicsConstants.MU_EARTH_KM3_S2 * (Units.KILOMETERS^3); 
    Re = PhysicsConstants.R_EARTH_KM * Units.KILOMETERS;     
    J2 = PhysicsConstants.J2_EARTH;
    J3 = PhysicsConstants.J3_EARTH;
    J4 = PhysicsConstants.J4_EARTH;
    
    % 2. Transform Position to ECEF
    % Vallado's r_I, r_J, r_K correspond to ECEF components for zonal harmonics
    r_ECEF = R_ECEF_from_ECI * r_ECI_meters;
    rI = r_ECEF(1);
    rJ = r_ECEF(2);
    rK = r_ECEF(3);
    r  = norm(r_ECEF);
    
    r2 = r*r;
    rK2 = rK*rK;
    
    % Initialize ECEF acceleration vector
    a_ECEF = [0; 0; 0];
    
    % 3. J2 Perturbation (Eq. 8-52)
    if J2_on
        coeffJ2 = -(3 * J2 * mu * Re^2) / (2 * r^5);
        
        a_ECEF = a_ECEF + coeffJ2 * [ ...
            rI * (1 - 5*rK2/r2); ...
            rJ * (1 - 5*rK2/r2); ...
            rK * (3 - 5*rK2/r2)  ...
        ];
    end
    
    % 4. J3 Perturbation (Eq. 8-53)
    if J3_on
        coeffJ3 = -(5 * J3 * mu * Re^3) / (2 * r^7);
        
        a_ECEF = a_ECEF + coeffJ3 * [ ...
            rI * (3*rK - 7*rK^3/r2); ...
            rJ * (3*rK - 7*rK^3/r2); ...
            rK * (6*rK2 - 7*rK^4/r2 - 0.6*r2) ...
        ];
    end
    
    % 5. J4 Perturbation (Eq. 8-54)
    if J4_on
        coeffJ4 = (15 * J4 * mu * Re^4) / (8 * r^7);
        
        a_ECEF = a_ECEF + coeffJ4 * [ ...
            rI * (1 - 14*rK2/r2 + 21*rK^4/r^4); ...
            rJ * (1 - 14*rK2/r2 + 21*rK^4/r^4); ...
            rK * (5 - 70/3*rK2/r2 + 21*rK^4/r^4) ...
        ];
    end
    
    % 6. Rotate Acceleration back to ECI
    accel_ECI_meters_s2 = R_ECEF_from_ECI' * a_ECEF;
end