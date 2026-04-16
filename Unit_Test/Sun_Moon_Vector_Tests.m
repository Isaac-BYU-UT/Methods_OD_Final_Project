classdef Sun_Moon_Vector_Tests < matlab.unittest.TestCase
    methods(Test)

        function test_Moon_Vector(testCase)
            [r_moon_rel_earth_ECI_meters_actual, ~, ~] = Forces.Vallado_moonPositionLowPrecision(2451545.0);

            [r_km, ~] = planetEphemeris(2451545.0, "Earth", "Moon");

            r_moon_rel_earth_ECI_meters_expected = r_km * Units.KILOMETERS;
            
            testCase.verifyEqual(r_moon_rel_earth_ECI_meters_actual(:), r_moon_rel_earth_ECI_meters_expected(:), 'AbsTol', 1e6);

        end

        function test_Sun_Vector(testCase)
            [r_sun_rel_earth_ECI_meters_actual, ~, ~] = Forces.Vallado_sunPositionLowPrecision(2451545.0);

            [r_km, ~] = planetEphemeris(2451545.0, "Earth", "Sun");

            r_sun_rel_earth_ECI_meters_expected = r_km * Units.KILOMETERS;
            
            testCase.verifyEqual(r_sun_rel_earth_ECI_meters_actual(:), r_sun_rel_earth_ECI_meters_expected(:), 'AbsTol', 1e7);

        end
    end
end