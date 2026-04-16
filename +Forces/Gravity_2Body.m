function accel_2Body_m_s2 = Gravity_2Body(r_ECI_meters)
    % GRAVITY_2BODY Calculates the acceleration due to a point mass Earth.
    %
    % Input:  r_ECI_meters     - 3x1 position vector [m]
    % Output: accel_2Body_m_s2 - 3x1 acceleration vector [m/s^2]
    
    r_mag_meters = norm(r_ECI_meters);
    mu = PhysicsConstants.MU_EARTH_KM3_S2 * (Units.KILOMETERS ^ 3);
    accel_2Body_m_s2 = -(mu / r_mag_meters^3) * r_ECI_meters;
end