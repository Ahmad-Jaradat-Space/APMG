function [grace_at_gps] = extractGRACEatGPS(u_vertical_grid, lat_grid, lon_grid, lat_gps, lon_gps)
% extractGRACEatGPS - Extract GRACE deformation at GPS station locations
%
% SYNTAX:
%   [grace_at_gps] = extractGRACEatGPS(u_vertical_grid, lat_grid, lon_grid, lat_gps, lon_gps)
%
% INPUT:
%   u_vertical_grid - GRACE vertical deformation grid [nlat x nlon] or [nlat x nlon x ntime]
%   lat_grid        - Latitude grid in degrees [-90 to 90] [nlat x nlon] or [nlat x 1]
%   lon_grid        - Longitude grid in degrees [-180 to 180] [nlat x nlon] or [1 x nlon]
%   lat_gps         - GPS station latitudes in degrees [nstations x 1]
%   lon_gps         - GPS station longitudes in degrees [nstations x 1]
%
% OUTPUT:
%   grace_at_gps    - GRACE deformation at GPS locations [nstations x ntime]
%
% METHOD:
%   - Uses bilinear interpolation for spatial interpolation
%   - Handles longitude wrapping (±180°)
%   - Validates GPS coordinates are within GRACE grid bounds
%   - Supports both single time and time series inputs
%
% NOTES:
%   - GRACE spatial resolution: ~300-500 km
%   - GPS measurements are point observations
%   - Some smoothing is expected due to resolution difference
%
% Author: GRACE Analysis Project
% Date: 2025

% Input validation
if nargin < 5
    error('All 5 input arguments are required');
end

% Handle different input formats for grids
if isvector(lat_grid) && isvector(lon_grid)
    [LON_GRID, LAT_GRID] = meshgrid(lon_grid, lat_grid);
elseif isequal(size(lat_grid), size(lon_grid))
    LAT_GRID = lat_grid;
    LON_GRID = lon_grid;
else
    error('lat_grid and lon_grid dimensions are inconsistent');
end

% Determine if we have time series or single time step
if ndims(u_vertical_grid) == 3
    [nlat, nlon, ntime] = size(u_vertical_grid);
    has_time_series = true;
else
    [nlat, nlon] = size(u_vertical_grid);
    ntime = 1;
    has_time_series = false;
end

% Validate grid dimensions
if ~isequal(size(LAT_GRID), [nlat, nlon]) || ~isequal(size(LON_GRID), [nlat, nlon])
    error('Grid dimensions do not match deformation grid');
end

% Validate GPS coordinates
nstations = length(lat_gps);
if length(lon_gps) ~= nstations
    error('lat_gps and lon_gps must have same length');
end

if any(lat_gps < -90 | lat_gps > 90)
    error('GPS latitudes must be between -90 and 90 degrees');
end

fprintf('Extracting GRACE deformation at %d GPS stations\n', nstations);
if has_time_series
    fprintf('Processing %d time steps\n', ntime);
end

%% Handle longitude wrapping
% Ensure longitudes are in consistent range
lon_gps_wrapped = mod(lon_gps + 180, 360) - 180; % [-180, 180]
LON_GRID_wrapped = mod(LON_GRID + 180, 360) - 180; % [-180, 180]

% Check if any GPS stations cross the date line
crosses_dateline = any(LON_GRID_wrapped(:) < -150) && any(LON_GRID_wrapped(:) > 150);
if crosses_dateline
    fprintf('Handling date line crossing...\n');
    % Convert to [0, 360] range to avoid interpolation issues
    LON_GRID_wrapped = mod(LON_GRID_wrapped, 360);
    lon_gps_wrapped = mod(lon_gps_wrapped, 360);
end

%% Validate GPS stations are within grid bounds
lat_min = min(LAT_GRID(:));
lat_max = max(LAT_GRID(:));
lon_min = min(LON_GRID_wrapped(:));
lon_max = max(LON_GRID_wrapped(:));

fprintf('Grid bounds: Lat [%.2f, %.2f], Lon [%.2f, %.2f]\n', ...
        lat_min, lat_max, lon_min, lon_max);

% Check if any GPS stations are outside bounds
outside_lat = lat_gps < lat_min | lat_gps > lat_max;
outside_lon = lon_gps_wrapped < lon_min | lon_gps_wrapped > lon_max;
outside_bounds = outside_lat | outside_lon;

if any(outside_bounds)
    warning('%d GPS stations are outside grid bounds:', sum(outside_bounds));
    for i = find(outside_bounds)'
        fprintf('  Station %d: Lat=%.3f, Lon=%.3f\n', i, lat_gps(i), lon_gps_wrapped(i));
    end
end

%% Initialize output array
grace_at_gps = zeros(nstations, ntime);

%% Extract deformation for each time step
for t = 1:ntime
    if has_time_series
        current_deformation = u_vertical_grid(:, :, t);
        if mod(t, max(1, floor(ntime/10))) == 0 || t == ntime
            fprintf('  Processing time step %d/%d (%.1f%%)\n', t, ntime, 100*t/ntime);
        end
    else
        current_deformation = u_vertical_grid;
    end
    
    % Check for valid data in current time step
    if all(all(isnan(current_deformation)))
        warning('All NaN values in deformation grid for time step %d', t);
        grace_at_gps(:, t) = NaN;
        continue;
    end
    
    % Extract values for each GPS station
    for s = 1:nstations
        try
            % Use MATLAB's interp2 for bilinear interpolation
            % Note: interp2 expects meshgrid format (X, Y, Z)
            interpolated_value = interp2(LON_GRID_wrapped, LAT_GRID, current_deformation, ...
                                       lon_gps_wrapped(s), lat_gps(s), 'linear');
            
            % Handle extrapolation for stations outside bounds
            if isnan(interpolated_value) && ~outside_bounds(s)
                % Try nearest neighbor if linear interpolation fails
                interpolated_value = interp2(LON_GRID_wrapped, LAT_GRID, current_deformation, ...
                                           lon_gps_wrapped(s), lat_gps(s), 'nearest');
            end
            
            grace_at_gps(s, t) = interpolated_value;
            
        catch ME
            warning('Interpolation failed for station %d at time %d: %s', s, t, ME.message);
            grace_at_gps(s, t) = NaN;
        end
    end
end

%% Quality assessment and statistics
fprintf('\nExtraction complete!\n');

% Count valid extractions
valid_extractions = ~isnan(grace_at_gps);
total_extractions = nstations * ntime;
n_valid = sum(valid_extractions(:));
success_rate = 100 * n_valid / total_extractions;

fprintf('Valid extractions: %d/%d (%.1f%%)\n', n_valid, total_extractions, success_rate);

% Per-station statistics
for s = 1:nstations
    station_valid = sum(valid_extractions(s, :));
    station_success = 100 * station_valid / ntime;
    
    if station_success < 50
        fprintf('Warning: Station %d has low success rate (%.1f%%)\n', s, station_success);
    end
end

% Overall statistics for valid data
if n_valid > 0
    valid_values = grace_at_gps(valid_extractions);
    
    fprintf('\nExtracted deformation statistics:\n');
    fprintf('  Mean: %.3f mm\n', mean(valid_values) * 1000);
    fprintf('  Std:  %.3f mm\n', std(valid_values) * 1000);
    fprintf('  Min:  %.3f mm\n', min(valid_values) * 1000);
    fprintf('  Max:  %.3f mm\n', max(valid_values) * 1000);
    fprintf('  RMS:  %.3f mm\n', sqrt(mean(valid_values.^2)) * 1000);
    
    % Seasonal analysis if time series
    if has_time_series && ntime >= 12
        fprintf('\nSeasonal analysis:\n');
        
        % Compute amplitude of seasonal signal for each station
        for s = 1:min(nstations, 5) % Show first 5 stations
            station_data = grace_at_gps(s, :);
            valid_idx = ~isnan(station_data);
            
            if sum(valid_idx) >= 12
                % Simple seasonal amplitude estimate
                seasonal_amp = (max(station_data(valid_idx)) - min(station_data(valid_idx))) / 2;
                fprintf('  Station %d seasonal amplitude: %.2f mm\n', s, seasonal_amp * 1000);
            end
        end
    end
    
else
    warning('No valid extractions - all values are NaN');
end

%% Additional validation
% Check for unrealistic values
if n_valid > 0
    max_abs_deformation = max(abs(valid_values)) * 1000; % in mm
    
    if max_abs_deformation > 200
        warning('Very large deformation values detected (max: %.1f mm)', max_abs_deformation);
    elseif max_abs_deformation < 0.01
        warning('Very small deformation values (max: %.3f mm)', max_abs_deformation);
    end
    
    % Check for consistent patterns across stations
    if nstations > 1 && has_time_series
        station_correlations = zeros(nstations, nstations);
        for i = 1:nstations
            for j = i+1:nstations
                valid_both = valid_extractions(i, :) & valid_extractions(j, :);
                if sum(valid_both) > 10
                    data_i = grace_at_gps(i, valid_both);
                    data_j = grace_at_gps(j, valid_both);
                    station_correlations(i, j) = corr(data_i', data_j');
                end
            end
        end
        
        % Report average inter-station correlation
        upper_tri = triu(station_correlations, 1);
        valid_corr = upper_tri(upper_tri ~= 0);
        if ~isempty(valid_corr)
            mean_correlation = mean(valid_corr);
            fprintf('Average inter-station correlation: %.3f\n', mean_correlation);
            
            if mean_correlation < 0.3
                warning('Low inter-station correlation - check for processing issues');
            end
        end
    end
end

fprintf('GPS coordinate range: Lat [%.3f, %.3f], Lon [%.3f, %.3f]\n', ...
        min(lat_gps), max(lat_gps), min(lon_gps), max(lon_gps));

end