function accel_2Body_km_s2 = Gravity_2Body(r_km)
    accel_2Body_km_s2 = gradient(Potential_2Body(r_km), r_km);
end

function U = Potential_2Body(r_km)
    r_mag_km = sqrt(r_km(1)^2 + r_km(2)^2 + r_km(3)^2);
    U = Constants.MU_EARTH_KM3_S2 / r_mag_km;
end