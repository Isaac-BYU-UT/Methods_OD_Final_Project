function [orthoDCM] = orthodcm_fast(DCM)
    % One iteration is usually enough for real-time drift correction
    % With some help from Gemini!
    orthoDCM = 0.5 * DCM * (3 * eye(3) - DCM' * DCM);
end