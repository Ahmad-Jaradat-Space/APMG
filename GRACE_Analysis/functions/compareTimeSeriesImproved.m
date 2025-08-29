function stats = compareTimeSeriesImproved(gps_daily, time_gps, grace_ts, grace_start_mjd, grace_end_mjd, time_grace)
%COMPARETIMESERIESIMPROVED Improved GPS-GRACE correlation with proper temporal alignment
%
% Based on peer-reviewed research for GPS-GRACE correlation:
% - Kurtenbach et al. (2012): Daily to monthly averaging techniques
% - Davis et al. (2004): GPS temporal correlation methods
% - Long et al. (2015): GRACE temporal filtering and alignment
%
% Input:
%   gps_daily - Daily GPS time series (meters)
%   time_gps - GPS time vector (MJD) 
%   grace_ts - GRACE monthly time series (meters)
%   grace_start_mjd - GRACE period start times (MJD)
%   grace_end_mjd - GRACE period end times (MJD)
%   time_grace - GRACE time vector (mid-points, MJD)
%
% Output:
%   stats - Structure with correlation statistics

% Initialize output structure
stats = struct();

% Apply low-pass filter to GPS data
fprintf('Applying low-pass filter to GPS data...\n');
gps_filtered = lowPassFilterGPS(gps_daily, time_gps);

% Average GPS over exact GRACE periods
fprintf('Averaging GPS to GRACE observation periods...\n');
[gps_monthly, valid_periods] = averageGPStoGRACEperiods(gps_filtered, time_gps, grace_start_mjd, grace_end_mjd);

% Extract corresponding GRACE values (only for valid periods)
gps_monthly_valid = gps_monthly(valid_periods);
grace_monthly = grace_ts(valid_periods);
time_common = time_grace(valid_periods);

% Remove any remaining NaN values
common_valid = ~isnan(gps_monthly_valid) & ~isnan(grace_monthly);
gps_common = gps_monthly_valid(common_valid);
grace_common = grace_monthly(common_valid);
time_final = time_common(common_valid);

n_common = length(gps_common);

if n_common < 5
    warning('Insufficient overlapping data points (%d) for reliable statistics', n_common);
    stats = struct('n_common', n_common, 'correlation', NaN, 'rmse', NaN, ...
                   'bias', NaN, 'nse', NaN, 'amplitude_gps', NaN, ...
                   'amplitude_grace', NaN, 'amplitude_ratio', NaN, 'best_lag_months', 0);
    return;
end

% Basic statistics
stats.n_common = n_common;
stats.bias = mean(grace_common - gps_common);
stats.rmse = sqrt(mean((grace_common - gps_common).^2));

% Correlation coefficient
if std(gps_common) > 1e-10 && std(grace_common) > 1e-10
    correlation_matrix = corrcoef(gps_common, grace_common);
    stats.correlation = correlation_matrix(1, 2);
else
    stats.correlation = NaN;
end

% Nash-Sutcliffe Efficiency
sse = sum((gps_common - grace_common).^2);
ssu = sum((gps_common - mean(gps_common)).^2);
stats.nse = 1 - sse/ssu;

% Amplitude comparison
stats.amplitude_gps = (max(gps_common) - min(gps_common)) / 2;
stats.amplitude_grace = (max(grace_common) - min(grace_common)) / 2;

if stats.amplitude_gps > 0
    stats.amplitude_ratio = stats.amplitude_grace / stats.amplitude_gps;
else
    stats.amplitude_ratio = NaN;
end

% Lag correlation analysis (up to 3 months)
if n_common >= 10
    max_lag = min(3, floor(n_common/3));  % Up to 3 months or 1/3 of data
    correlations = zeros(2*max_lag + 1, 1);
    lags = -max_lag:max_lag;
    
    for i = 1:length(lags)
        lag = lags(i);
        if lag == 0
            correlations(i) = stats.correlation;
        elseif lag > 0
            % GRACE leads GPS
            if length(grace_common) > lag && length(gps_common) > lag
                r = corrcoef(grace_common(1:end-lag), gps_common(1+lag:end));
                correlations(i) = r(1,2);
            end
        else
            % GPS leads GRACE
            lag_pos = -lag;
            if length(gps_common) > lag_pos && length(grace_common) > lag_pos
                r = corrcoef(gps_common(1:end-lag_pos), grace_common(1+lag_pos:end));
                correlations(i) = r(1,2);
            end
        end
    end
    
    [max_corr, best_idx] = max(abs(correlations));
    stats.best_lag_months = lags(best_idx);
    stats.max_correlation = correlations(best_idx);
else
    stats.best_lag_months = 0;
    stats.max_correlation = stats.correlation;
end

% Quality assessment
fprintf('GPS-GRACE Temporal Alignment Results:\n');
fprintf('  Common periods: %d\n', n_common);
fprintf('  Correlation: %.3f\n', stats.correlation);
fprintf('  Best lag correlation: %.3f (lag: %+d months)\n', stats.max_correlation, stats.best_lag_months);
fprintf('  RMSE: %.3f mm\n', stats.rmse * 1000);
fprintf('  Amplitude ratio: %.3f\n', stats.amplitude_ratio);

end