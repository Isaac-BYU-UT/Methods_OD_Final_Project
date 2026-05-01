function EGM96_Data = get_2020_EGM96_Data()

    persistent EGM96_Data_cached

    if isempty(EGM96_Data_cached)
    raw_egm96_coefficients = importdata("EOP_data/egm96_to360_ascii.txt");
    egm96_table = array2table(raw_egm96_coefficients,"VariableNames",{'l','m','Cbar_l_m','Sbar_l_m','std_C_n_m','std_S_n_m'}); % Note, the Original EGM data uses 'n' instead of 'l'.
    
    table_20_20_EGM96 = egm96_table(1:228,:);

    %% Constants
    L_max = 20; % For 20x20 EGM 96.
    
    % Pi_Normalization_Matrix = Tools.getNormalizationWeights(L_max);
    
    C_coeff_matrix = NaN(L_max+1, L_max+1);
    S_coeff_matrix = NaN(L_max+1, L_max+1);
    for i = 1:height(table_20_20_EGM96)
        L = table_20_20_EGM96.l(i) + 1; % Note we always have to add 1 for MATLAB!
        m = table_20_20_EGM96.m(i) + 1; % Note we always have to add 1 for MATLAB!
        C_coeff_matrix(L,m) = table_20_20_EGM96.Cbar_l_m(i); % / Pi_Normalization_Matrix(L,m);
        S_coeff_matrix(L,m) = table_20_20_EGM96.Sbar_l_m(i); % / Pi_Normalization_Matrix(L,m);
    end

    EGM96_Data_cached.C_coeff_matrix = C_coeff_matrix;
    EGM96_Data_cached.S_coeff_matrix = S_coeff_matrix;
    % EGM96_Data_cached.Pi_Normalization_Matrix = Pi_Normalization_Matrix;
    EGM96_Data_cached.L_max = L_max;

    end
        EGM96_Data = EGM96_Data_cached;

    
    
end