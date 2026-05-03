for i_run_case = 1:6
    EKF.run_ekf(i_run_case,false)
end

EKF.run_ekf_short_arc;

% EKF.run_ekf(6,true)


%% Plot Batches
% 
% file_names = {'Results/5_2_2026_All_Cases/sorensen.mat', ...
%               'Results/5_2_2026_All_Cases_B/sorensen.mat',...
%               'Results/5_2_2026_All_Cases_C/sorensen.mat',...
%               'Results/5_2_2026_More_Crosstrack_Noise/sorensen.mat',...
%               'Results/5_2_2026_Even_More_Y/sorensen.mat' };


file_names = {'Results/5_2_2026_All_Cases/sorensen.mat', ...
              'Results/5_2_2026_All_Cases_B/sorensen.mat',...
              'Results/5_2_2026_All_Cases_C/sorensen.mat'}

Visuals.plot3DPositions_multi(file_names)
