function Pi_matrix = getNormalizationWeights(L_max)
    % Generated with help from Gemini
    % Generates a matrix of Pi_lm coefficients based on Eq. 8-22
    % Pi_lm = sqrt( (l+m)! / ( (l-m)! * k * (2l+1) ) )
    
    Pi_matrix = NaN(L_max + 1, L_max + 1);
    
    for L = 0:L_max
        for m = 0:L
            % Define k based on m
            if m == 0
                k = 1;
            else
                k = 2;
            end
            
            % Compute using factorials as per Eq. 8-22
            num = factorial(L + m);
            den = factorial(L - m) * k * (2*L + 1);
            
            Pi_matrix(L+1, m+1) = sqrt(num / den);
        end
    end
end