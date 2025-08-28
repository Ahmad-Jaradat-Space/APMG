function [A] = get_designmatrx_harmonic(t, l, p, omega, q, t_jump, t_quake)
% GET_DESIGNMATRX_HARMONIC - Constructs the design matrix for polynomial and harmonic approximation
%
% INPUTS:
%   t     - time vector (Nx1), e.g., decimal years
%   l     - observations (Nx1), not used in design matrix but kept for interface compatibility
%   p     - degree of polynomial
%   omega - base angular frequency (e.g. 2*pi for annual terms)
%   q     - vector of frequency multipliers (e.g. [1 2] for annual + semiannual)
% t_jump: one (or more) jump dates in MJD
% t_quake_list: vector of earthquake dates (in MJD)

% OUTPUT:

%   A     - design matrix (NxM), where M = p+1 + 2*length(q)



% Normalize time
t_jump = t_jump - t(1);
t_quake = t_quake - t(1);
t = t-t(1);
t_norm = t / max(t);
n = length(t);

% initial design matrix
A = ones(n, 1); 


% polynomial terms:
A_poly = t_norm.^[0 : p];


% harmonic terms: cos(q*omega*t), sin(q*omega*t)
%for i = 1:length(q)
%    A = [A, cos(q(i)*omega*t), sin(q(i)*omega*t)];
%end

% harmonic terms:
if isempty(q)
    A_har = [];
else
    A_har = [cos(q.*omega.*t), sin(q.*omega.*t)];
end

% jump column:
% J = zeros(size(t));
% J = [];
% if nargin > 5 && ~isempty(t_jump)
%     J = double(t >= reshape(t_jump, 1, []));
%     idx = J(end,:) == 0;
%     J(:,idx) = [];
% end

% Indicator columns for each jump/quake
J = zeros(length(t), length(t_jump));
for k = 1:length(t_jump)
    if k < length(t_jump)
        J(:,k) = (t >= t_jump(k)) & (t < t_jump(k+1));
    else
        J(:,k) = (t >= t_jump(k));
    end
end


% Log colomn:
Log = [];
if nargin > 6 && ~isempty(t_quake)
    delta = t - reshape(t_quake, 1, []);
    mask = delta >= 0;
    Log = zeros(size(delta));
    Log(mask) = log(1 + delta(mask));
    idx = mask(end,:) == 0;
    Log(:,idx) = [];
end



% final design matrix:
A = [A_poly, A_har, J, Log];



%if nargin > 5 && ~isempty(t_jump)
%    J = double(t + t(1) >= t_jump); % t + t(1) = original t
%    A = [A_poly, A_har, J];
%else
%    A = [A_poly, A_har];
%end

end