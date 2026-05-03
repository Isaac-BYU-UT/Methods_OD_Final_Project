function IC_Sat_Epoch = loadSatInitialConditions_Short_Arc()
    IC_Sat_Epoch.position_ECI_meters = [ -2.141405286772074e+06; -6.782912806990056e+06; -2.495347836732329e+05]; % m
    IC_Sat_Epoch.velocity_ECI_meters_per_second = [7.155970076085581e+03; -2.252009394465575e+03; 34.521991629986779]; % m/s
    IC_Sat_Epoch.epoch_date_time_UTC = datetime(2018, 3, 23, 8, 55, 3.0, 'TimeZone','UTC'); % Don't add seconds here, add it later on.
end


% Updated with EKF_Results_F_Final_6D_Batched_2026-05-02_16-34-02

% I've Arrived at 5 Days, my state is: 
%    1.0e+06 *
% 
%   -2.141405294214197
%   -6.782912722768536
%   -0.249535385625627
%    0.007155970095932
%   -0.002252009394539
%    0.000034521856807

   % Based on results from
   % EKF_Results_F_Final_6D_Batched_2026-05-02_14-44-56.mat

   % RIGHT BEFORE!!! At index: 2134 of ekf.t_obs (432360 seconds), the one
   % right below it!



   