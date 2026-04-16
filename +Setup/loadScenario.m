function Scenario = loadScenario(case_name, days_to_fit)
Scenario.time_start_seconds = 0; % seconds since epoch  
Scenario.dt_seconds = 60; % seconds

switch days_to_fit

    case '1 Hour'
        Scenario.time_end_seconds = 3600; % seconds

    case '6 Hours'
        Scenario.time_end_seconds = 21600; % seconds

    case '1 Day'
        Scenario.time_end_seconds = 86400; % seconds

    case '3 Days'
        Scenario.time_end_seconds = 3 * 86400; % seconds

    case '6 Days'
        Scenario.time_end_seconds = 6 * 86400; % seconds
end

Scenario.tspan_seconds = Scenario.time_start_seconds:Scenario.dt_seconds:Scenario.time_end_seconds;

Scenario.range_on = true;
Scenario.range_rate_on = true;
Scenario.Atoll_on = true;
Scenario.Diego_Garcia_on = true;
Scenario.Arecibo_on = true;

switch case_name
    case 'A' % Range Only, All Sensors
        Scenario.range_rate_on = false;
        return;

    case 'B' % Range-Rate Only, All Sensors
        Scenario.range_on = false;
        return;

    case 'C' % Atoll Only, All Data Types
        Scenario.Diego_Garcia_on = false;
        Scenario.Arecibo_on = false;
        return;

    case 'D' % Diego Garcia Only, All Data Types
        Scenario.Atoll_on = false;
        Scenario.Arecibo_on = false;
        return;

    case 'E' % Arecibo Only, All Data Types
        Scenario.Atoll_on = false;
        Scenario.Diego_Garcia_on = false;
        return;

    case 'F' % Fit the Long Arc, All Data, All Stations
        return;

    case 'G' % Fit the Short Arc, Only The Last Day of Data for All Sensors
        Scenario.tspan_seconds = (Scenario.time_end_seconds - 86400):Scenario.dt_seconds:Scenario.time_end_seconds;
        return;

end