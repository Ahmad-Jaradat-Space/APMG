function [h] = globalTest(vTPv, n, u, sigma0)
% Short description of the function.
%  globalTest performs a global model test (Chi-square test) for LSA.
% Input:
%   vTPv - residuals' quadratic form (v'Pv)
%   n    - number of observations
%   u    - number of unknowns (i.e., parameters estimated)
% Output:
%   h    - test result: 1 = passed (model is valid), 0 = rejected
 

% Degrees of freedom
f = n - u;

% Confidence level
alpha = 0.05;

% Critical value bounds
chi2_lower = chi2inv(alpha/2, f);
chi2_upper = chi2inv(1 - alpha/2, f);

% Test statistic
T = vTPv / sigma0^2;

% Decision
h = (T >= chi2_lower) && (T <= chi2_upper);


end

