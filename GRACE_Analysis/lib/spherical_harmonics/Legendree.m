function [P]= Legendree (lat,L)
% Legendree - Associated Legendre polynomials (batch computation)
%
% SYNTAX:
%   [P] = Legendree(lat, L)
%
% INPUT:
%   lat - Colatitude in degrees [scalar or vector]
%   L   - Maximum degree
%
% OUTPUT:
%   P   - Associated Legendre polynomials [L+1 x L+1]
%         P(n+1, m+1) contains P_nm for degree n, order m
%
% METHOD:
%   Efficient computation of normalized associated Legendre functions
%   following geodetic conventions. Uses recursive formulation for
%   computational efficiency.
%
% NOTES:
%   - Input latitude is in degrees (consistent with main_fun_new.m)
%   - Outputs normalized Legendre functions for spherical harmonics
%   - More efficient than individual pnm() calls for batch computation
%
% REFERENCE:
%   Standard geodetic recursive formulation for Associated Legendre polynomials
%
% Author: Original Satellite Geodesy Assignment
% Adapted: GRACE Analysis Project, 2025

% Convert latitude to trigonometric functions
c = cosd(lat);  % cosine of colatitude
s = sind(lat);  % sine of colatitude

% Initialize Legendre polynomial matrix
P = zeros(L+1, L+1);

% Recursive computation of Associated Legendre polynomials
for n = 0:1:L
    for m = 0:1:n
        
        % Normalization factors
        an = sqrt((2*n + 1) / (2*n));
        bn = sqrt(2*n + 1);
        cn = sqrt((2*n + 1) / ((n - m) * (n + m)));
        dn = sqrt(2*n - 1);
        en = sqrt(((n - m - 1) * (n + m - 1)) / (2*n - 3));
        
        % Recursive computation based on degree and order
        if (n == 0 && m == 0)
            % Base case: P_00 = 1
            P(n+1, m+1) = 1;
            
        elseif (n == 1 && m == 1)
            % Special case: P_11 = sqrt(3) * sin(colatitude)
            P(n+1, m+1) = sqrt(3) * s;
            
        elseif (n == m)
            % Diagonal terms: P_nn
            P(n+1, m+1) = an * s * P(n, m+1);
            
        elseif (n - m == 1)
            % Off-diagonal terms: P_(m+1,m)
            P(n+1, m+1) = bn * c * P(n, m+1);
            
        else
            % General recursive formula
            P(n+1, m+1) = cn * (dn * c * P(n, m+1) - en * P(n-1, m+1));
        end
    end
end

end