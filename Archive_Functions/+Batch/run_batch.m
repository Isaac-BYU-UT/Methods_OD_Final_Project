function run_batch()
% Batch Least Squares Orbital Determination
% Main driver script following EKF module conventions
% Processes all measurements at once using batch least squares algorithm

    clear; clc;

    %% SETUP
    i_case = 6;
    case_names = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

    i_scenario = 2;
    scenario_config = {'Accel', 'Final_1D', 'Final_3D', 'Final_6D', 'HW5', '24Dynamics'};

    [S, ENV] = Setup.loadSettings(case_names{i_case}, scenario_config{i_scenario}, false, false);

    %% INITIALIZE BATCH STRUCTURE
    batch.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
    batch.time_struct_epoch = ENV.time_struct_epoch;
    batch.print_updates = true;
    batch.debug_on = false;
    batch.ode_type = 'ode113';
    batch.f_updates = 999;

    [S, batch] = Batch.initialize_batch(S, ENV, batch);

    if batch.print_updates
        disp('Initial State:'); disp(transpose(batch.X0_nominal));
        disp('A priori Covariance:'); disp(batch.P_cov_0);
    end

    %% ITERATIVE BATCH LEAST SQUARES (UNTIL CONVERGENCE)
    max_iterations = 10;
    convergence_tol = 1e-6;  % Convergence tolerance on norm of state correction
    iteration = 0;
    converged = false;
    dX_batch = zeros(batch.N_states, 1); % Initialize vector for state update
    P_cov_batch = batch.P_cov_0;
    
    W_apriori = P_cov_batch \ eye(batch.N_states);
    Information_Matrix = W_apriori;
    Residual_Vector = 0;

    while ~converged && iteration < max_iterations
        iteration = iteration + 1;

        if batch.print_updates
            fprintf('\n=== ITERATION %d ===\n', iteration);
        end


        %% PROPAGATE ALL STATES AND STATE TRANSITION MATRICES
        if batch.print_updates
            disp('Propagating states to all observation times...');
        end

        [X_states_all, STM_all, y_ode_all] = Batch.propagate_states_batch(batch, S, ENV);

        if batch.print_updates
            fprintf('Propagation complete: %d states computed\n', batch.N_obs);
        end

        %% BUILD DESIGN MATRIX (H MATRIX) FOR ALL MEASUREMENTS
        if batch.print_updates
            disp('Building design matrix for all measurements...');
        end

        [H_batch, curr_meas_all] = Batch.build_measurement_matrix(batch, S, ENV, X_states_all, STM_all);

        if batch.print_updates
            fprintf('Design matrix built: %d x %d\n', size(H_batch, 1), size(H_batch, 2));
        end

        %% BUILD OBSERVATION VECTOR WITH LIGHT-TIME CORRECTIONS
        if batch.print_updates
            disp('Building observation vector with light-time corrections...');
        end

        [y_batch, curr_meas_all] = Batch.build_observation_vector(batch, S, ENV, X_states_all, curr_meas_all);

        if batch.print_updates
            fprintf('Observation vector built: %d measurements\n', length(y_batch));
        end

        %% SOLVE BATCH LEAST SQUARES
        if batch.print_updates
            disp('Solving batch least squares normal equations...');
        end

        [dX_batch, P_cov_batch, batch] = Batch.solve_batch_least_squares(batch, S, H_batch, y_batch, curr_meas_all);

        if batch.print_updates
            disp('Solution obtained.');
        end

        %% CHECK CONVERGENCE
        dX_norm = norm(dX_batch);
        if batch.print_updates
            fprintf('Norm of state correction: %e\n', dX_norm);
            fprintf('Convergence tolerance: %e\n', convergence_tol);
        end

        if dX_norm < convergence_tol
            converged = true;
            if batch.print_updates
                fprintf('CONVERGED after %d iterations!\n\n', iteration);
            end
        end

        %% UPDATE STATE ESTIMATE FOR NEXT ITERATION
        batch.X0_nominal = batch.X0_nominal + dX_batch;
        batch.X_input = batch.X0_nominal;

        if batch.print_updates
            fprintf('Updated state estimate. New state:\n');
            disp(transpose(batch.X0_nominal));
        end

    end

    if ~converged && batch.print_updates
        fprintf('WARNING: Did not converge after %d iterations\n', max_iterations);
    end

   

end
