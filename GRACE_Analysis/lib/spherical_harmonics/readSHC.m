function [cnm, snm] = readSHC(file, nmax)
% readSHC - Read spherical harmonic coefficients with enforced nmax
% 
% INPUTS:
%   file - coefficient data [n, m, cnm, snm]
%   nmax - maximum degree (optional, defaults to max degree in file)
%
% OUTPUTS:
%   cnm, snm - coefficient matrices of size (nmax+1, nmax+1)

if nargin < 2
    nmax = max(file(:, 1));
end

n = file(:, 1);
m = file(:, 2);
cnm = zeros(nmax + 1, nmax + 1);
snm = zeros(nmax + 1, nmax + 1);

for i = 1:size(file, 1)
    degree = n(i) + 1;
    order = m(i) + 1;
    % Only include coefficients up to nmax
    if degree <= nmax + 1 && order <= nmax + 1
        cnm(degree, order) = file(i, 3);
        snm(degree, order) = file(i, 4);
    end
end
end