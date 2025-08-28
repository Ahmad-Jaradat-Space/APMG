% read the SHCs
function [cnm, snm] = readSHC(file)
% readSHC - Read spherical harmonic coefficients from GRACE data matrix
%
% INPUT:
%   file - Matrix with columns [n, m, cnm, snm] where n=degree, m=order
%
% OUTPUT:
%   cnm - Cosine coefficients matrix (n+1 x m+1)
%   snm - Sine coefficients matrix (n+1 x m+1)

if isempty(file)
    error('Input file matrix is empty');
end

if size(file, 2) < 4
    error('Input matrix must have at least 4 columns [n, m, cnm, snm]');
end

% Extract degree and order
n = file(:, 1);
m = file(:, 2);

% Find maximum degree to initialize matrices
nmax = max(n);

% Initialize coefficient matrices with zeros
cnm = zeros(nmax + 1, nmax + 1);
snm = zeros(nmax + 1, nmax + 1);

% Fill coefficient matrices
for i = 1:size(file, 1)
    degree = n(i) + 1;  % Convert to MATLAB 1-based indexing
    order = m(i) + 1;   % Convert to MATLAB 1-based indexing
    
    % Bounds checking
    if degree > 0 && order > 0 && degree <= nmax + 1 && order <= degree
        cnm(degree, order) = file(i, 3);
        snm(degree, order) = file(i, 4);
    end
end

end