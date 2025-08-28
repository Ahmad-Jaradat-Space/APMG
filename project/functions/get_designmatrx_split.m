function [A_com] = get_designmatrx_split(t, l, p, omega, q, t_jump, t_quake, t_split)
% GET_DESIGNMATRX_SPLIT - Construct a block-diagonal design matrix that fits
% two independent models before and after a split time (e.g., an earthquake).
%
% INPUTS:
%   t        - time vector (Nx1), e.g., in MJD (Modified Julian Date)
%   l        - observation vector (Nx1), e.g., east/north/up component
%   p        - polynomial degree
%   omega    - base angular frequency for harmonics
%   q        - vector of frequency multipliers (e.g., [1 2] for annual, semiannual)
%   t_jump   - vector of jump times (MJD)
%   t_quake  - vector of earthquake times (MJD)
%   t_split  - time at which to split the model (MJD)
%
% OUTPUT:
%   A_com    - Block-diagonal design matrix combining models before/after split

    % Find indices for data before and after the split
    idx_1 = t < t_split;    % Before split
    idx_2 = t >= t_split;   % After split

    % Construct design matrix for the "before" segment
    A_1 = get_designmatrx_harmonic(t(idx_1), l(idx_1), p, omega, q, t_jump, t_quake);

    % Construct design matrix for the "after" segment
    A_2 = get_designmatrx_harmonic(t(idx_2), l(idx_2), p, omega, q, t_jump, t_quake);

    % Build the block-diagonal matrix:
    % [ A_1   0  ]
    % [  0   A_2 ]
    % Rows are for all times; columns are parameters for before and after.
    A_com = [A_1, zeros(size(A_1,1),size(A_2,2));    
             zeros(size(A_2,1),size(A_1,2)), A_2];   
end
