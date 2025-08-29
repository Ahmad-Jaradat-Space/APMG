function [h_n, l_n, k_n] = loadLoveNumbers(nmax, model)
% Try to load PREM load Love numbers from MAT file; fallback to legacy
try
    mat_path = fullfile('data','love_numbers','PREM_Love_Numbers_n60.mat');
    if exist(mat_path, 'file')
        data = load(mat_path);
        % Attempt common variable names; adjust if needed
        candidates_h = {'h','h_n','h_load','h_load_n'};
        candidates_k = {'k','k_n','k_load','k_load_n'};
        candidates_l = {'l','l_n','l_load','l_load_n'};
        h_all = [];
        k_all = [];
        l_all = [];
        for fn = candidates_h
            if isfield(data, fn{1}), h_all = data.(fn{1}); break; end
        end
        for fn = candidates_k
            if isfield(data, fn{1}), k_all = data.(fn{1}); break; end
        end
        for fn = candidates_l
            if isfield(data, fn{1}), l_all = data.(fn{1}); break; end
        end
        if isempty(h_all) || isempty(k_all) || isempty(l_all)
            error('Love number arrays not found in MAT file');
        end
        % Ensure column vectors
        h_all = h_all(:); k_all = k_all(:); l_all = l_all(:);
        N = min([nmax+1, length(h_all), length(k_all), length(l_all)]);
        h_n = h_all(1:N);
        k_n = k_all(1:N);
        l_n = l_all(1:N);
        % Pad if needed
        if length(h_n) < nmax+1
            h_n(end+1:nmax+1) = h_n(end);
            k_n(end+1:nmax+1) = k_n(end);
            l_n(end+1:nmax+1) = l_n(end);
        end
        return;
    end
catch
    % Fall through to legacy values
end

% Fallback legacy values (not recommended scientifically)
h_n = zeros(nmax + 1, 1);
l_n = zeros(nmax + 1, 1);
k_n = zeros(nmax + 1, 1);
h_n(1) = 0;
l_n(1) = 0;
k_n(1) = 0;
if nmax >= 1
    h_n(2) = 0; l_n(2) = 0; k_n(2) = 0.021;
end
for n = 2:nmax
    if n == 2
        h_n(n+1) = 0.6149;
        l_n(n+1) = 0.0839;
        k_n(n+1) = 0.3020;
    else
        h_n(n+1) = 0.6149 * (1 + 5/(2*n+1) - 7/(2*n+1)^2);
        l_n(n+1) = 0.0839 * (1 + 3/(2*n+1) - 5/(2*n+1)^2);
        k_n(n+1) = 0.3020 * (1 + 2/(2*n+1) - 3/(2*n+1)^2);
    end
end
end