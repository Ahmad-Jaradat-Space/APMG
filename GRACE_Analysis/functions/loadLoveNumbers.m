function [h_n, l_n, k_n] = loadLoveNumbers(nmax, model)
% Load PREM load Love numbers (h', l', k') to degree nmax, consistent with
% fully normalized GRACE coefficients. If the MAT contains tidal Love numbers
% (h, l, k), convert to load using 1 + k_tidal - h_tidal = 1 + k_load ⇒ k' = k - h.

mat_path = fullfile('data','love_numbers','PREM_Love_Numbers_n60.mat');
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


end