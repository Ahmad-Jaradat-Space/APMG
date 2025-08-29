function [gps_monthly, grace_periods] = averageGPStoGRACEperiods(gps_daily, time_gps, grace_start_mjd, grace_end_mjd)
%AVERAGEGPSTOGRACEPERIODS Average daily GPS data over exact GRACE observation periods
%
% Based on peer-reviewed research (Kurtenbach et al. 2012, Davis et al. 2004):
% "As a reference method for combining daily data into monthly series, 
% researchers use the arithmetic mean of the days exactly covered by 
% each GRACE monthly solution"
%
% Input:
%   gps_daily - Daily GPS time series (meters)
%   time_gps - GPS time vector (MJD)
%   grace_start_mjd - Start times of GRACE observation periods (MJD)
%   grace_end_mjd - End times of GRACE observation periods (MJD)
%
% Output:
%   gps_monthly - GPS data averaged over GRACE periods
%   grace_periods - Logical array indicating which periods have valid data

n_grace_periods = length(grace_start_mjd);
gps_monthly = NaN(n_grace_periods, 1);
grace_periods = false(n_grace_periods, 1);

% Ensure column vectors
gps_daily = gps_daily(:);
time_gps = time_gps(:);

for i = 1:n_grace_periods
    % Find GPS observations within this GRACE period
    period_mask = time_gps >= grace_start_mjd(i) & time_gps <= grace_end_mjd(i);
    
    if sum(period_mask) >= 3  % Need at least 3 days for reliable average
        gps_in_period = gps_daily(period_mask);
        
        % Remove outliers using robust statistics (2.5 sigma)
        median_val = median(gps_in_period, 'omitnan');
        mad_val = mad(gps_in_period, 1);  % Mean absolute deviation
        outlier_threshold = 2.5 * mad_val;
        
        clean_mask = abs(gps_in_period - median_val) <= outlier_threshold;
        
        if sum(clean_mask) >= 1
            gps_monthly(i) = mean(gps_in_period(clean_mask), 'omitnan');
            grace_periods(i) = true;
        end
    end
end


end