function [stats] = compareTimeSeries(gps_ts, grace_ts, time_gps, time_grace)
% compareTimeSeries - Statistical comparison of GPS and GRACE time series
%
% SYNTAX:
%   [stats] = compareTimeSeries(gps_ts, grace_ts, time_gps, time_grace)
%
% INPUT:
%   gps_ts     - GPS vertical displacement time series [n_gps x 1]
%   grace_ts   - GRACE vertical deformation time series [n_grace x 1]  
%   time_gps   - GPS time stamps in MJD [n_gps x 1]
%   time_grace - GRACE time stamps in MJD [n_grace x 1]
%
% OUTPUT:
%   stats      - Structure containing comparison statistics:
%              .correlation      - Pearson correlation coefficient
%              .correlation_p    - P-value for correlation significance
%              .rmse            - Root Mean Square Error [m]
%              .nse             - Nash-Sutcliffe Efficiency
%              .bias            - Mean bias (GRACE - GPS) [m]
%              .mae             - Mean Absolute Error [m]
%              .nrmse           - Normalized RMSE (%)
%              .amplitude_gps   - GPS seasonal amplitude [m]
%              .amplitude_grace - GRACE seasonal amplitude [m]
%              .amplitude_ratio - GRACE/GPS amplitude ratio
%              .phase_lag       - Phase lag between signals [days]
%              .n_common        - Number of common observations
%              .time_common     - Common time vector [MJD]
%              .gps_interp      - Interpolated GPS values
%              .grace_interp    - Interpolated GRACE values
%
% METHOD:
%   1. Time series alignment using linear interpolation
%   2. Statistical metrics computation using existing functions
%   3. Seasonal analysis using FFT-based methods
%   4. Uncertainty estimation and significance testing
%
% REFERENCES:
%   - Fu & Freymueller (2012) for GPS-GRACE comparison methodology
%   - Nash & Sutcliffe (1970) for NSE metric
%
% Author: GRACE Analysis Project
% Date: 2025

% Add paths to existing statistical functions
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));

% Input validation
if nargin < 4
    error('All 4 input arguments are required');
end

% Ensure all inputs are column vectors
gps_ts = gps_ts(:);
grace_ts = grace_ts(:);
time_gps = time_gps(:);
time_grace = time_grace(:);

% Check dimensions
if length(gps_ts) ~= length(time_gps)
    error('GPS time series and time vector must have same length');
end

if length(grace_ts) ~= length(time_grace)
    error('GRACE time series and time vector must have same length');
end

fprintf('Comparing GPS and GRACE time series\n');
fprintf('GPS data points: %d\n', length(gps_ts));
fprintf('GRACE data points: %d\n', length(grace_ts));

%% Step 1: Remove NaN values from individual series
gps_valid = ~isnan(gps_ts) & ~isnan(time_gps);
grace_valid = ~isnan(grace_ts) & ~isnan(time_grace);

gps_ts_clean = gps_ts(gps_valid);
time_gps_clean = time_gps(gps_valid);
grace_ts_clean = grace_ts(grace_valid);
time_grace_clean = time_grace(grace_valid);

fprintf('After NaN removal: GPS=%d, GRACE=%d\n', length(gps_ts_clean), length(grace_ts_clean));

if length(gps_ts_clean) < 2 || length(grace_ts_clean) < 2
    error('Insufficient data points after NaN removal');
end

%% Step 2: Determine common time period
time_start = max(min(time_gps_clean), min(time_grace_clean));
time_end = min(max(time_gps_clean), max(time_grace_clean));

fprintf('Common time period: MJD %.1f to %.1f\n', time_start, time_end);

if time_end <= time_start
    error('No overlapping time period between GPS and GRACE data');
end

% Create common time grid (monthly sampling)
time_common = time_start:30.44:time_end; % ~30.44 days per month
time_common = time_common(:);

fprintf('Common time grid: %d points\n', length(time_common));

%% Step 3: Interpolate both time series to common grid
try
    % Remove duplicate time points in GPS data
    [time_gps_unique, gps_unique_idx] = unique(time_gps_clean, 'stable');
    gps_ts_unique = gps_ts_clean(gps_unique_idx);
    
    % Remove duplicate time points in GRACE data
    [time_grace_unique, grace_unique_idx] = unique(time_grace_clean, 'stable');
    grace_ts_unique = grace_ts_clean(grace_unique_idx);
    
    % Check for sufficient unique points
    if length(time_gps_unique) < 2 || length(time_grace_unique) < 2
        error('Insufficient unique time points for interpolation (GPS=%d, GRACE=%d)', ...
              length(time_gps_unique), length(time_grace_unique));
    end
    
    % GPS interpolation with unique time points
    gps_interp = interp1(time_gps_unique, gps_ts_unique, time_common, 'linear', NaN);
    
    % GRACE interpolation with unique time points
    grace_interp = interp1(time_grace_unique, grace_ts_unique, time_common, 'linear', NaN);
    
catch ME
    error('Interpolation failed: %s', ME.message);
end

%% Step 4: Find common valid observations
common_valid = ~isnan(gps_interp) & ~isnan(grace_interp);
n_common = sum(common_valid);

fprintf('Common valid observations: %d\n', n_common);

if n_common < 4
    error('Insufficient common observations (%d) for statistical analysis', n_common);
elseif n_common < 6
    warning('Limited common observations (%d) - results may have high uncertainty', n_common);
end

% Extract common data
gps_common = gps_interp(common_valid);
grace_common = grace_interp(common_valid);
time_common_valid = time_common(common_valid);

%% Step 5: Basic statistical metrics
fprintf('\nComputing statistical metrics...\n');

% Initialize stats structure
stats = struct();
stats.n_common = n_common;
stats.time_common = time_common_valid;
stats.gps_interp = gps_common;
stats.grace_interp = grace_common;

% Mean bias
stats.bias = mean(grace_common - gps_common);

% Mean Absolute Error
stats.mae = mean(abs(grace_common - gps_common));

% Root Mean Square Error using existing rms.m function
try
    stats.rmse = rms(grace_common - gps_common);
catch
    % Fallback if rms function not available
    stats.rmse = sqrt(mean((grace_common - gps_common).^2));
end

% Normalized RMSE (percentage of GPS standard deviation)
gps_std = std(gps_common);
if gps_std > 0
    stats.nrmse = 100 * stats.rmse / gps_std;
else
    stats.nrmse = Inf;
end

% Correlation calculation - using corrcoef which is always available
correlation_matrix = corrcoef(gps_common, grace_common);
stats.correlation = correlation_matrix(1, 2);

% Simple significance test (rough approximation)
if abs(stats.correlation) < 1
    t_stat = stats.correlation * sqrt((n_common - 2) / (1 - stats.correlation^2));
    stats.correlation_p = 2 * (1 - tcdf(abs(t_stat), n_common - 2));
else
    stats.correlation_p = 0; % Perfect correlation
end

% Nash-Sutcliffe Efficiency using existing NSE.m function
try
    stats.nse = NSE(gps_common, grace_common);
catch
    % Fallback calculation
    ss_res = sum((gps_common - grace_common).^2);
    ss_tot = sum((gps_common - mean(gps_common)).^2);
    stats.nse = 1 - (ss_res / ss_tot);
end

fprintf('Basic statistics computed:\n');
fprintf('  Correlation: %.3f (p=%.4f)\n', stats.correlation, stats.correlation_p);
fprintf('  RMSE: %.3f mm\n', stats.rmse * 1000);
fprintf('  NSE: %.3f\n', stats.nse);
fprintf('  Bias: %.3f mm\n', stats.bias * 1000);

%% Step 6: Seasonal analysis
fprintf('\nPerforming seasonal analysis...\n');

% Estimate seasonal amplitude for each time series (requires at least 6 points for meaningful analysis)
if n_common >= 6
    [stats.amplitude_gps, phase_gps] = estimateSeasonalAmplitude(gps_common, time_common_valid);
    [stats.amplitude_grace, phase_grace] = estimateSeasonalAmplitude(grace_common, time_common_valid);
else
    % Fallback for limited data
    stats.amplitude_gps = (max(gps_common) - min(gps_common)) / 2;
    stats.amplitude_grace = (max(grace_common) - min(grace_common)) / 2;
    phase_gps = 0;
    phase_grace = 0;
end

% Amplitude ratio
if stats.amplitude_gps > 0
    stats.amplitude_ratio = stats.amplitude_grace / stats.amplitude_gps;
else
    stats.amplitude_ratio = NaN;
end

% Phase lag (in days)
stats.phase_lag = (phase_grace - phase_gps) * 365.25 / (2 * pi);

% Ensure phase lag is in [-182.5, 182.5] days
if stats.phase_lag > 182.5
    stats.phase_lag = stats.phase_lag - 365.25;
elseif stats.phase_lag < -182.5
    stats.phase_lag = stats.phase_lag + 365.25;
end

fprintf('Seasonal analysis results:\n');
fprintf('  GPS amplitude: %.2f mm\n', stats.amplitude_gps * 1000);
fprintf('  GRACE amplitude: %.2f mm\n', stats.amplitude_grace * 1000);
fprintf('  Amplitude ratio: %.2f\n', stats.amplitude_ratio);
fprintf('  Phase lag: %.1f days\n', stats.phase_lag);

%% Step 7: Quality assessment
fprintf('\nQuality assessment:\n');

% Correlation significance
if stats.correlation_p < 0.001
    fprintf('  Correlation: Highly significant (p < 0.001)\n');
elseif stats.correlation_p < 0.01
    fprintf('  Correlation: Very significant (p < 0.01)\n');
elseif stats.correlation_p < 0.05
    fprintf('  Correlation: Significant (p < 0.05)\n');
else
    fprintf('  Correlation: Not significant (p = %.3f)\n', stats.correlation_p);
end

% NSE interpretation
if stats.nse > 0.75
    fprintf('  NSE: Excellent model performance (%.3f)\n', stats.nse);
elseif stats.nse > 0.50
    fprintf('  NSE: Good model performance (%.3f)\n', stats.nse);
elseif stats.nse > 0.30
    fprintf('  NSE: Acceptable model performance (%.3f)\n', stats.nse);
else
    fprintf('  NSE: Poor model performance (%.3f)\n', stats.nse);
end

% RMSE in context
seasonal_range_gps = max(gps_common) - min(gps_common);
rmse_percentage = 100 * stats.rmse / seasonal_range_gps;
fprintf('  RMSE: %.2f mm (%.1f%% of GPS seasonal range)\n', ...
        stats.rmse * 1000, rmse_percentage);

% Amplitude ratio assessment
if abs(stats.amplitude_ratio - 1) < 0.2
    fprintf('  Amplitude agreement: Good (ratio = %.2f)\n', stats.amplitude_ratio);
elseif abs(stats.amplitude_ratio - 1) < 0.5
    fprintf('  Amplitude agreement: Fair (ratio = %.2f)\n', stats.amplitude_ratio);
else
    fprintf('  Amplitude agreement: Poor (ratio = %.2f)\n', stats.amplitude_ratio);
end

end

%% Helper function for seasonal amplitude estimation
function [amplitude, phase] = estimateSeasonalAmplitude(data, time_mjd)
% Estimate seasonal amplitude and phase using least squares fitting
%
% INPUT:
%   data     - Time series data
%   time_mjd - Time in MJD
%
% OUTPUT:
%   amplitude - Seasonal amplitude
%   phase     - Phase in radians

    % Remove mean
    data_demean = data - mean(data);
    
    % Convert time to decimal years
    time_year = 1858.0 + time_mjd / 365.25; % Approximate conversion
    
    % Design matrix for annual and semi-annual terms
    t = time_year - time_year(1); % Time relative to start
    
    A = [ones(length(t), 1), t, ...
         cos(2*pi*t), sin(2*pi*t), ...
         cos(4*pi*t), sin(4*pi*t)];
    
    try
        % Least squares solution
        params = A \ data_demean;
        
        % Annual components
        cos_coeff = params(3);
        sin_coeff = params(4);
        
        % Amplitude and phase
        amplitude = sqrt(cos_coeff^2 + sin_coeff^2);
        phase = atan2(sin_coeff, cos_coeff);
        
    catch
        % Fallback: simple range-based amplitude
        amplitude = (max(data) - min(data)) / 2;
        phase = 0;
    end
end