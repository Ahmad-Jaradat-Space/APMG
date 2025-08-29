function [h_n, l_n, k_n] = loadLoveNumbers(nmax, model)
% Load PREM load Love numbers (h', l', k') to degree nmax, consistent with
% fully normalized GRACE coefficients. If the MAT contains tidal Love numbers
% (h, l, k), convert to load using 1 + k_tidal - h_tidal = 1 + k_load ⇒ k' = k - h.

mat_path = fullfile('data','love_numbers','PREM_Love_Numbers_n60.mat');
if ~exist(mat_path, 'file')
    error('Love numbers MAT not found: %s. Please provide PREM load Love numbers.', mat_path);
end

S = load(mat_path);
% Try common variable names
candidates_h = {'h_load','h_load_n','h_prime','h_n','h'};
candidates_k = {'k_load','k_load_n','k_prime','k_n','k'};
candidates_l = {'l_load','l_load_n','l_prime','l_n','l'};

h_all = [];
k_all = [];
l_all = [];
for fn = candidates_h
    if isfield(S, fn{1}), h_all = S.(fn{1}); break; end
end
for fn = candidates_k
    if isfield(S, fn{1}), k_all = S.(fn{1}); break; end
end
for fn = candidates_l
    if isfield(S, fn{1}), l_all = S.(fn{1}); break; end
end
if isempty(h_all) || isempty(k_all) || isempty(l_all)
    error('Could not find h/k/l arrays in %s', mat_path);
end

% Column vectors and truncate/pad to nmax+1
h_all = h_all(:); k_all = k_all(:); l_all = l_all(:);
N = min([nmax+1, length(h_all), length(k_all), length(l_all)]);
h_n = h_all(1:N);
k_n = k_all(1:N);
l_n = l_all(1:N);
if length(h_n) < nmax+1
    h_n(end+1:nmax+1) = h_n(end);
    k_n(end+1:nmax+1) = k_n(end);
    l_n(end+1:nmax+1) = l_n(end);
end

% Detect whether k_n is tidal (typically positive around degree 2) and convert to load
% For load Love numbers, k'_n (n>=2) is typically negative (~ -0.3 for n=2).
idx = 3:min(nmax+1, 6); % degrees 2..5
if any(k_n(idx) > 0)
    % Assume tidal k; convert to load: k' = k - h
    k_n = k_n - h_n;
end

% Basic sanity check on degree-2
deg2 = 3; % index for n=2
if deg2 <= length(h_n)
    % Expect h2' ~ 0.5-0.7, k2' ~ -0.45..-0.1, so 1+k2' ~ 0.3..0.9
    if ~(h_n(deg2) > 0.4 && h_n(deg2) < 0.8)
        warning('Suspicious h''_2 = %.3f; verify MAT content.', h_n(deg2));
    end
    if ~(k_n(deg2) < 0 && (1 + k_n(deg2)) > 0.2 && (1 + k_n(deg2)) < 1.2)
        warning('Suspicious k''_2 = %.3f; transformed to load? Verify inputs.', k_n(deg2));
    end
end

end