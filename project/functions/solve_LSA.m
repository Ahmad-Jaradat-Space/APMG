function [out] = solve_LSA(A, l, sigma0)
% SOLVE_LSA - Performs Least Squares Adjustment
%
% INPUTS:
%   A       - Design matrix (NxM)
%   l       - Observation vector (Nx1)
%   sigma0  - A-priori standard deviation of unit weight (scalar)
%
% OUTPUT:
%   out - structure with:
%    out.x_est     - Estimated parameters
%    out.Sigma_X   - Covariance matrix of estimated parameters
%    out.l_hat     - Adjusted observations
%    out.Sigma_L   - Covariance matrix of adjusted observations
%    out.v         - Residuals
%    out.sigma0hat - Estimated variance factor (global test)

% Number of observations and unknowns
n = size(A, 1);
m = size(A, 2);

% Covariance matrix of observations
Sigma_L = sigma0^2 * eye(n);

% Normal matrix and right-hand side
N = A' * A;
u = A' * l;

% Solve normal equations
x_est = N \ u;

% Compute residuals
l_hat = A * x_est;
v = l - l_hat;

% Variance factor estimate
sigma0hat2 = (v' * v) / (n - m);  % Global model test: estimated variance of unit weight

% Covariance of estimated parameters
Sigma_X = sigma0hat2 * inv(N);

% Covariance of adjusted observations
Sigma_L_adj = A * Sigma_X * A';

% Output
out.x_est = x_est;
out.Sigma_X = Sigma_X;
out.l_hat = l_hat;
out.Sigma_L = Sigma_L_adj;
out.v = v;
out.sigma0hat = sigma0hat2;

end
