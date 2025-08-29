function debug_p056_break()
close all; clc;
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'time_utils'));

stn = 'P056';
file = fullfile('data','gps',[stn '.tenv3']);
if ~exist(file, 'file')
    fprintf('File not found: %s\n', file);
    return;
end

S = load_tenv3(file);
if ~isfield(S,'t') || ~isfield(S,'up') || isempty(S.t) || isempty(S.up)
    fprintf('Missing data fields in %s\n', file);
    return;
end

t_mjd = S.t(:);
y = S.up(:);

valid = isfinite(t_mjd) & isfinite(y);
t_mjd = t_mjd(valid);
y = y(valid);

if numel(y) < 100
    fprintf('Insufficient data points: %d\n', numel(y));
    return;
end

% Convert time to decimal years
try
    t_year = mjd2decyear(t_mjd);
catch
    t_year = 1858 + t_mjd/365.25; % fallback
end

% Heuristic unit check: convert to mm if likely meters
yrange = prctile(abs(y),[5 95]);
if yrange(2) < 1 % likely meters
    y_mm = y*1000;
    ulabel = 'mm';
else
    y_mm = y;
    ulabel = 'mm';
end

N = numel(y_mm);

% Single linear fit (baseline)
A0 = [ones(N,1) t_year];
coeff0 = A0\y_mm;
res0 = y_mm - A0*coeff0;
SSE0 = sum(res0.^2);

% Two-segment fit with one break (scan)
min_seg = max(30, round(0.1*N));
idx_min = min_seg;
idx_max = N - min_seg;
SSE_best = inf;
best_k = NaN;
coeff1_best = [NaN;NaN];
coeff2_best = [NaN;NaN];

for k = idx_min:idx_max
    t1 = t_year(1:k); y1 = y_mm(1:k);
    t2 = t_year(k+1:end); y2 = y_mm(k+1:end);
    A1 = [ones(numel(t1),1) t1];
    A2 = [ones(numel(t2),1) t2];
    c1 = A1\y1; c2 = A2\y2;
    SSE = sum((y1 - A1*c1).^2) + sum((y2 - A2*c2).^2);
    if SSE < SSE_best
        SSE_best = SSE;
        best_k = k;
        coeff1_best = c1;
        coeff2_best = c2;
    end
end

% Chow F-test
p_restr = 2;           % parameters in restricted (single) model
p_unrestr = 4;         % parameters in unrestricted (two segments)
num = (SSE0 - SSE_best)/(p_unrestr - p_restr);
den = SSE_best/(N - p_unrestr);
F = num/den;
pval = 1 - fcdf(F, p_unrestr - p_restr, N - p_unrestr);

% Report
break_time = t_year(best_k);
break_date = datestr(datenum(0,1,1) + (break_time - floor(break_time))*365.25 + datenum(floor(break_time),0,0)); %#ok<DATNM>

slope1 = coeff1_best(2); % mm/year
slope2 = coeff2_best(2); % mm/year

fprintf('Station %s structural break analysis\n', stn);
fprintf('  N = %d points\n', N);
fprintf('  Best break index: %d / %d\n', best_k, N);
fprintf('  Best break time (decimal year): %.4f\n', break_time);
fprintf('  Slope before break:  %.3f %s/yr\n', slope1, ulabel);
fprintf('  Slope after break:   %.3f %s/yr\n', slope2, ulabel);
fprintf('  SSE single: %.3e, SSE two-segment: %.3e\n', SSE0, SSE_best);
fprintf('  Chow F(%d,%d) = %.2f, p = %.3g\n', p_unrestr - p_restr, N - p_unrestr, F, pval);

if pval < 0.05
    fprintf('Conclusion: Significant break detected (p < 0.05). Two linear trends with an offset are appropriate.\n');
else
    fprintf('Conclusion: No statistically significant break at 5%%; single linear may suffice.\n');
end

end 