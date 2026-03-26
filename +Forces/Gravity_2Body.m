function accel_2Body_km_s2 = Gravity_2Body(r_km)
    % r_km must be in ECI coordinates and in km for this function to work properly. The output acceleration will be in km/s^2.
    r_mag_km = sqrt(r_km(1)^2 + r_km(2)^2 + r_km(3)^2);
    U = Constants.MU_EARTH_KM3_S2 / r_mag_km;
    accel_2Body_km_s2 = gradient(U, r_km);
end
