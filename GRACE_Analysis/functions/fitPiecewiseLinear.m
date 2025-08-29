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

% Create detrended series and fitted values
y_detrended = y;
y_fitted = y;
y_detrended(1:break_idx) = y1 - (slope1 * x1 + intercept1);
y_detrended(break_idx:end) = y2 - (slope2 * x2 + intercept2);
y_fitted(1:break_idx) = slope1 * x1 + intercept1;
y_fitted(break_idx:end) = slope2 * x2 + intercept2;

result.v = y_detrended;
result.lest = y_fitted;
result.slope1 = slope1;
result.slope2 = slope2;
result.intercept1 = intercept1;
result.intercept2 = intercept2;
result.break_point = break_point;
result.break_idx = break_idx;
end 