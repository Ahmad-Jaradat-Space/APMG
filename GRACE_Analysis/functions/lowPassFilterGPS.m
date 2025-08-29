function gps_filtered = lowPassFilterGPS(gps_ts, time_gps)
%LOWPASSFILTERGPS Apply low-pass filter to GPS time series before averaging
%
% Based on peer-reviewed research (Long et al. 2015):
% "Apply a low-pass Butterworth filter with 12 cpy cut-off frequency 
% to remove any remaining signal with sub-monthly periods from daily 
% time series"
%
% Input:
%   gps_ts - GPS time series (daily, meters)
%   time_gps - GPS time vector (MJD)
%
% Output:
%   gps_filtered - Low-pass filtered GPS time series

% Check for sufficient data
if length(gps_ts) < 50
    warning('GPS time series too short for reliable filtering');
    gps_filtered = gps_ts;
    return;
end

% Remove NaN values and interpolate gaps < 7 days
valid_mask = ~isnan(gps_ts);
if sum(valid_mask) < 0.8 * length(gps_ts)
    warning('GPS time series has too many gaps (>20%%)');
end

% Interpolate short gaps
gps_clean = gps_ts;
gap_start = [];
gap_end = [];
in_gap = false;

for i = 1:length(valid_mask)
    if ~valid_mask(i) && ~in_gap
        gap_start = [gap_start, i];
        in_gap = true;
    elseif valid_mask(i) && in_gap
        gap_end = [gap_end, i-1];
        in_gap = false;
    end
end

% Close final gap if needed
if in_gap
    gap_end = [gap_end, length(valid_mask)];
end

% Fill short gaps (<=7 days) with linear interpolation
for i = 1:length(gap_start)
    gap_length = gap_end(i) - gap_start(i) + 1;
    if gap_length <= 7  % Max 7 days
        if gap_start(i) > 1 && gap_end(i) < length(gps_ts)
            % Linear interpolation
            x1 = gap_start(i) - 1;
            x2 = gap_end(i) + 1;
            y1 = gps_clean(x1);
            y2 = gps_clean(x2);
            
            for j = gap_start(i):gap_end(i)
                weight = (j - x1) / (x2 - x1);
                gps_clean(j) = y1 + weight * (y2 - y1);
            end
        end
    end
end

% Calculate sampling frequency from median time difference
dt_median = median(diff(time_gps), 'omitnan');  % days
fs = 1 / dt_median;  % samples per day

% Design Butterworth filter
% 12 cycles per year = 12/365.25 cycles per day
fc = 12 / 365.25;  % cutoff frequency (cycles per day)
fn = fc / (fs/2);  % normalized frequency

% Ensure normalized frequency is valid
if fn >= 1
    warning('Cutoff frequency too high relative to sampling rate');
    fn = 0.9;  % Use 90% of Nyquist
end

try
    [b, a] = butter(4, fn, 'low');  % 4th order Butterworth
    
    % Apply zero-phase filtering to preserve temporal relationships
    gps_filtered = filtfilt(b, a, gps_clean);
    
    
catch ME
    warning('Filter design failed: %s. Using unfiltered data.', ME.message);
    gps_filtered = gps_clean;
end

end