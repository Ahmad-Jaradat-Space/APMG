function [out] = fitPolynomial(t, y, p, sigma0)
% fitPolynomial fits a polynomial of degree p to coordinate time series using LSA.
%
% Input:
%   t       - time vector 
%   y       - coordinate vector (east, north, up)
%   p       - degree of the polynomial
%   sigma0  - standard deviation of observations (mean)
%
% Output:
%   out.xest        - estimated parameter vector
%   out.lest        - adjusted observations
%   out.Cxx         - covariance matrix of the parameters
%   out.Cll         - covariance matrix of adjusted observations
%   out.v           - residuals
%   out.global_test - result of global model test

% number of observations
nObs = length(t);
t_norm = t / max(t);
% design matrix
A = zeros(nObs, p+1);
for i = 0:p
    A(:, i+1) = t_norm.^i;
end

%A = t_norm.^[0 : p];

% Observation vector
l = y;

% Covariance matrix of observations 
Cll = sigma0^2 * eye(nObs);
P = inv(Cll);

% Normal matrix and right-hand side
N = A' * P * A;
rhs = A' * P * l;

% Solve normal equations
out.xest = N \ rhs;


% Adjusted observations
out.lest = A * out.xest;

% residual vector
out.v = l - out.lest;

% Covariance matrices
vTPv = out.v' * P * out.v;
out.s0_hat = sqrt(vTPv / (nObs - (p + 1)));
out.Cxx = out.s0_hat^2 * inv(N);
out.Cll = A * out.Cxx * A';



end