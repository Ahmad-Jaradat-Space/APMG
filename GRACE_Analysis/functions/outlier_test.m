function [outliers, tau, r, Qv_diag] = outlier_test(A, v, sigma_value, isGlobalTestAccepted, sigLevel)
% Performs Baarda outlier test or Pope's test depending on global model test result
%
% Inputs:
%   A - design matrix (n x m)
%   v - residuals (n x 1)
%   sigma_value - sigma0 if global test passed, sqrt(sigma0hat) if failed
%   isGlobalTestAccepted - boolean: true = Baarda test, false = Pope test
%   sigLevel - significance level (e.g., 0.05 or 0.01)
%
% Outputs:
%   outliers - logical vector of detected outliers
%   tau      - test statistic for each residual
%   r        - redundancy numbers (same as diag(Qv))
%   Qv_diag  - alias of r for clarity

    n = length(v);
    m = size(A, 2);
    Qv = eye(n) - A * ((A' * A) \ A');
    r = diag(Qv);  % redundancy numbers

    % Test statistic (standardized residuals)
    tau = abs(v) ./ (sigma_value * sqrt(r));

    % Compute critical value
    if isGlobalTestAccepted
        p = 1 - sigLevel / 2;
        threshold = sqrt(2) * erfinv(2 * p - 1);
    else
        % Pope test: t-distribution with (n - m - 1) DOF
        threshold = tinv(1 - sigLevel / 2, n - m - 1);
    end

    % Identify outliers
    outliers = tau > threshold;

    % For completeness
    Qv_diag = r;
end
