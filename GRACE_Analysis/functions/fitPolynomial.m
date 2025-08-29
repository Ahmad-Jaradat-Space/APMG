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

function result = fitPiecewiseLinear(x, y, break_point)
% Fit two linear segments with a break point
% x: time vector (MJD or decimal years)
% y: data vector
% break_point: break time (same units as x)
% Returns: result with v (detrended), slopes, intercepts, break_point

% Find break index
[~, break_idx] = min(abs(x - break_point));

% Split data
x1 = x(1:break_idx);
y1 = y(1:break_idx);
x2 = x(break_idx:end);
y2 = y(break_idx:end);

% Fit first segment
p1 = polyfit(x1, y1, 1);
slope1 = p1(1);
intercept1 = p1(2);

% Fit second segment  
p2 = polyfit(x2, y2, 1);
slope2 = p2(1);
intercept2 = p2(2);

% Create detrended series
y_detrended = y;
y_detrended(1:break_idx) = y1 - (slope1 * x1 + intercept1);
y_detrended(break_idx:end) = y2 - (slope2 * x2 + intercept2);

result.v = y_detrended;
result.slope1 = slope1;
result.slope2 = slope2;
result.intercept1 = intercept1;
result.intercept2 = intercept2;
result.break_point = break_point;
result.break_idx = break_idx;
end