function [H_batch, curr_meas_all] = build_measurement_matrix(batch, S, ENV, X_states_all, STM_all)
% Build design matrix (H) for all measurements in batch least squares
% Following the same measurement model as EKF

    N_states = batch.N_states;
    N_obs = batch.N_obs;
    t_obs = batch.t_obs;

    % Pre-allocate design matrix (2*N_obs rows for range and range-rate)
    H_batch = zeros(2 * N_obs, N_states);

    % Pre-allocate measurement structure array
    curr_meas_all = struct('station_id', zeros(N_obs, 1), ...
                           'y_obs_meters', zeros(2, N_obs), ...
                           'R', cell(N_obs, 1), ...
                           'r_stn_ECI_m', zeros(3, N_obs), ...
                           'v_stn_ECI_m_s', zeros(3, N_obs), ...
                           'H_matrix', zeros(2, N_states, N_obs));

    % Build H matrix row by row for each measurement
    for i = 1:N_obs
        % Get current measurement info
        station_id = S.ref_data.Actual_Measurements.station_id(i);
        
        % Store station ID
        curr_meas_all(i).station_id = station_id;

        % Get observation (apparent range and range rate)
        curr_meas_all(i).y_obs_meters = [S.ref_data.Actual_Measurements.apparent_range_km(i); ...
                                         S.ref_data.Actual_Measurements.apparent_range_rate_km_s(i)] ...
                                        * Units.KILOMETERS;

        % Get measurement covariance
        curr_meas_all(i).R = S.Stations(station_id).Covariance;

        % Interpolate EOP at measurement time
        date_time_at_meas = S.IC_Sat_Epoch.epoch_date_time_UTC + seconds(t_obs(i));
        EOP_at_meas = Tools.interpolate_EOP(date_time_at_meas, ENV.EOP_IERS, ENV.EOP_Celestrak);

        % Get station state in ECI at measurement time
        [curr_meas_all(i).r_stn_ECI_m, curr_meas_all(i).v_stn_ECI_m_s] = ...
            Tools.ECEF_ECI_Converter( ...
                S.Stations(station_id).position_ECEF_meters, ...
                zeros(3, 1), date_time_at_meas, "ECEF_to_ECI", EOP_at_meas);

        % Get propagated state at this time
        r_sat_ECI_m = X_states_all(1:3, i);
        v_sat_ECI_m_s = X_states_all(4:6, i);

        % Compute H matrix (partial derivatives) for this measurement
        H_i = Measurements.Compute_H_matrix(r_sat_ECI_m, v_sat_ECI_m_s, ...
                                            curr_meas_all(i).r_stn_ECI_m, curr_meas_all(i).v_stn_ECI_m_s);

        % Store H matrix for this measurement
        curr_meas_all(i).H_matrix = H_i;

        % Fill in design matrix rows for this measurement
        H_batch(2*i-1:2*i, :) = H_i * STM_all(:,:,i);
    end

end
