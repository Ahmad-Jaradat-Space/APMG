function [is_significant, t_values, crit_value] = significance_test(x_est, Cov_X, jump_idx, alpha, globalTestAccepted)
% significance_test
% Significance test for jump parameters in LSA results.
% 
% INPUTS:
%   x_est             : Estimated parameter vector (from LSA)
%   Cov_X             : Covariance matrix of x_est (from LSA, e.g., inv(A'*A) * sigma0^2)
%   jump_idx          : Indices of jump parameters in x_est (vector)
%   alpha             : Significance level 
%   globalTestAccepted: Boolean, true if global model test was accepted (use normal distribution),
%                                      false if declined (use Student's t, empirical sigma)
%
% OUTPUTS:
%   is_significant : Logical vector, true if jump is significant
%   t_values       : Test statistic for each jump
%   crit_value     : Critical value used (z_alpha or t_alpha)

if nargin < 5
    globalTestAccepted = true; % Default: assume global test passed
end

n_jumps = length(jump_idx);
is_significant = false(n_jumps,1);
t_values = zeros(n_jumps,1);

if globalTestAccepted
    crit_value = norminv(1-alpha/2,0,1); % ≈ 1.96 for alpha=0.05
else
    % Use Student's t for unknown variance, degrees of freedom: length(x_est)-length(jump_idx)
    dof = max(1, length(x_est) - length(jump_idx));
    crit_value = tinv(1-alpha/2, dof);
end

for j = 1:n_jumps
    idx = jump_idx(j);
    xk = x_est(idx);
    sigma_xk = sqrt(abs(Cov_X(idx, idx))); % abs for numerical safety
    
    t_values(j) = abs(xk / sigma_xk);
    is_significant(j) = t_values(j) >= crit_value;
end

end
