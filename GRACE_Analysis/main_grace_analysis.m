function main_grace_analysis()
% main_grace_analysis - Complete GRACE-GPS analysis workflow
%
% PURPOSE:
%   Comprehensive analysis of elastic crustal deformation using GRACE 
%   satellite data and GPS observations following Wahr et al. (1998) 
%   and Fu & Freymueller (2012) methodologies.
%
% WORKFLOW:
%   1. Load Love numbers (PREM model)
%   2. Process GRACE files with C20/degree-1 corrections
%   3. Load and process GPS time series data
%   4. Convert GRACE to vertical deformation grids
%   5. Extract GRACE at GPS station locations
%   6. Compare GPS and GRACE time series
%   7. Generate comprehensive results and visualizations
%
% REFERENCES:
%   - Wahr et al. (1998) JGR, hydrological deformation theory
%   - Fu & Freymueller (2012) JGR, GPS-GRACE comparison methodology
%   - Farrell (1972) Rev Geophys, elastic loading theory
%
% Author: Hamza Jaradat
% Date: 2025

%% Initialize analysis
close all; clear; clc;

fprintf('==========================================\n');
fprintf('GRACE-GPS Crustal Deformation Analysis\n');
fprintf('Based on Wahr et al. (1998) & Fu & Freymueller (2012)\n');
fprintf('==========================================\n\n');

% Record analysis start time
analysis_start = tic;

% Add all necessary paths
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib', 'gps_processing'));

%% Configuration parameters
config = struct();

% GRACE processing parameters
config.nmax = 60;                          % Maximum SH degree for GRACE
config.earth_model = 'PREM';               % Earth model for Love numbers
config.grace_dir = 'data/grace';           % Updated path for current structure
config.c20_file = 'data/aux/C20_RL05.txt'; % Updated path for current structure
config.deg1_file = 'data/aux/deg1_coef.txt'; % Updated path for current structure

% GPS processing parameters
config.gps_data_dir = 'data/gps';          % Updated path for current structure  
config.gps_coords_file = 'data/gps/GPSLatLong.tenv3'; % Updated to use actual coordinate file
config.detrend_method = 'linear';          % Linear or polynomial detrending

% Spatial analysis parameters
config.lat_range = [30, 50];               % Analysis latitude range [deg]
config.lon_range = [-130, -110];           % Analysis longitude range [deg]
config.grid_res = 1.0;                     % Grid resolution [degrees]

% Output configuration
config.save_intermediate = true;           % Save intermediate results
config.generate_plots = true;              % Create visualizations
config.output_dir = 'output';

% Create output directory
if ~exist(config.output_dir, 'dir')
    mkdir(config.output_dir);
end

%% Step 1: Load Love numbers
fprintf('Step 1: Loading Love numbers (%s model)...\n', config.earth_model);
step1_start = tic;

try
    [h_n, l_n, k_n] = loadLoveNumbers(config.nmax, config.earth_model);
    
    fprintf('  Love numbers loaded successfully\n');
    fprintf('  h_2 = %.5f, l_2 = %.5f, k_2 = %.5f\n', h_n(3), l_n(3), k_n(3));
    fprintf('  Step 1 completed in %.2f seconds\n\n', toc(step1_start));
    
catch ME
    error('Step 1 failed: %s', ME.message);
end

%% Step 2: Process GRACE time series
fprintf('Step 2: Processing GRACE coefficient files...\n');
step2_start = tic;

try
    [cnm_ts, snm_ts, time_mjd_grace] = processGRACEfiles(...
        config.grace_dir, config.c20_file, config.deg1_file);
    
    n_months = length(time_mjd_grace);
    fprintf('  Processed %d monthly GRACE solutions\n', n_months);
    fprintf('  Time range: MJD %.1f to %.1f\n', min(time_mjd_grace), max(time_mjd_grace));
    fprintf('  Step 2 completed in %.2f seconds\n\n', toc(step2_start));
    
    % Save intermediate results
    if config.save_intermediate
        save(fullfile(config.output_dir, 'grace_coefficients.mat'), ...
             'cnm_ts', 'snm_ts', 'time_mjd_grace', '-v7.3');
    end
    
catch ME
    error('Step 2 failed: %s', ME.message);
end

%% Step 3: Load and process GPS data
fprintf('Step 3: Loading GPS time series data...\n');
step3_start = tic;

try
    % Load GPS coordinates
    gps_coords = load(config.gps_coords_file);
    station_names = gps_coords(:, 1);  % Assuming first column is station ID
    lat_gps = gps_coords(:, 2);        % Latitude
    lon_gps = gps_coords(:, 3);        % Longitude
    
    n_stations = length(lat_gps);
    fprintf('  Found %d GPS stations\n', n_stations);
    
    % Initialize GPS data storage
    gps_data = cell(n_stations, 1);
    time_mjd_gps = cell(n_stations, 1);
    
    % Process each GPS station
    valid_stations = false(n_stations, 1);
    
    for i = 1:n_stations
        try
            % Construct filename (assuming .tenv3 format)
            station_file = fullfile(config.gps_data_dir, ...
                                  sprintf('station_%04d.tenv3', station_names(i)));
            
            if exist(station_file, 'file')
                % Load time series using existing function
                [gps_ts] = load_tenv3(station_file);
                
                if ~isempty(gps_ts) && size(gps_ts, 1) > 50  % Minimum 50 observations
                    % Extract vertical component and time
                    time_mjd_gps{i} = gps_ts(:, 1);  % MJD time
                    elevation = gps_ts(:, 4);        % Vertical component
                    
                    % Detrend GPS time series
                    if strcmp(config.detrend_method, 'linear')
                        elevation_detrended = fitPolynomial(time_mjd_gps{i}, elevation, 1);
                    else
                        elevation_detrended = elevation;  % No detrending
                    end
                    
                    gps_data{i} = elevation_detrended;
                    valid_stations(i) = true;
                    
                    if mod(i, 10) == 0 || i == n_stations
                        fprintf('    Processed %d/%d stations\n', i, n_stations);
                    end
                end
            end
            
        catch
            % Skip problematic stations
            valid_stations(i) = false;
        end
    end
    
    % Filter to valid stations only
    lat_gps = lat_gps(valid_stations);
    lon_gps = lon_gps(valid_stations);
    gps_data = gps_data(valid_stations);
    time_mjd_gps = time_mjd_gps(valid_stations);
    n_valid_stations = sum(valid_stations);
    
    fprintf('  Successfully processed %d valid GPS stations\n', n_valid_stations);
    fprintf('  Step 3 completed in %.2f seconds\n\n', toc(step3_start));
    
catch ME
    error('Step 3 failed: %s', ME.message);
end

%% Step 4: Create analysis grid and convert GRACE to deformation
fprintf('Step 4: Converting GRACE to vertical deformation...\n');
step4_start = tic;

try
    % Create analysis grid
    lat_vec = config.lat_range(1):config.grid_res:config.lat_range(2);
    lon_vec = config.lon_range(1):config.grid_res:config.lon_range(2);
    [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
    
    % Convert to radians for spherical harmonics
    theta_grid = (90 - lat_grid) * pi / 180;  % Colatitude in radians
    lambda_grid = lon_grid * pi / 180;        % Longitude in radians
    
    fprintf('  Analysis grid: %d x %d points\n', length(lat_vec), length(lon_vec));
    fprintf('  Geographic bounds: %.1f°-%.1f°N, %.1f°-%.1f°W\n', ...
            config.lat_range, -config.lon_range(2), -config.lon_range(1));
    
    % Initialize deformation time series
    [nlat, nlon] = size(lat_grid);
    u_vertical_ts = zeros(nlat, nlon, n_months);
    
    % Process each time step
    fprintf('  Computing vertical deformation for %d time steps:\n', n_months);
    
    for t = 1:n_months
        % Extract coefficients for this time step
        cnm = cnm_ts(:, :, t);
        snm = snm_ts(:, :, t);
        
        % Convert to vertical deformation using Wahr et al. (1998)
        u_vertical = graceToVerticalDeformation(cnm, snm, theta_grid, lambda_grid, h_n, k_n);
        
        u_vertical_ts(:, :, t) = u_vertical;
        
        if mod(t, 12) == 0 || t == n_months
            fprintf('    Processed %d/%d time steps (%.1f%%)\n', t, n_months, 100*t/n_months);
        end
    end
    
    fprintf('  Step 4 completed in %.2f seconds\n\n', toc(step4_start));
    
    % Save intermediate results
    if config.save_intermediate
        save(fullfile(config.output_dir, 'grace_deformation.mat'), ...
             'u_vertical_ts', 'lat_grid', 'lon_grid', 'time_mjd_grace', '-v7.3');
    end
    
catch ME
    error('Step 4 failed: %s', ME.message);
end

%% Step 5: Extract GRACE at GPS locations
fprintf('Step 5: Extracting GRACE deformation at GPS stations...\n');
step5_start = tic;

try
    % Extract GRACE values at GPS locations for all time steps
    grace_at_gps_ts = extractGRACEatGPS(u_vertical_ts, lat_grid, lon_grid, lat_gps, lon_gps);
    
    fprintf('  Extracted GRACE values for %d stations and %d time steps\n', ...
            n_valid_stations, n_months);
    fprintf('  Step 5 completed in %.2f seconds\n\n', toc(step5_start));
    
catch ME
    error('Step 5 failed: %s', ME.message);
end

%% Step 6: Compare GPS and GRACE time series
fprintf('Step 6: Statistical comparison of GPS and GRACE time series...\n');
step6_start = tic;

try
    % Initialize results storage
    comparison_stats = cell(n_valid_stations, 1);
    
    fprintf('  Comparing time series for %d stations:\n', n_valid_stations);
    
    for i = 1:n_valid_stations
        % Get GPS and GRACE data for this station
        gps_ts = gps_data{i};
        grace_ts = grace_at_gps_ts(i, :)';  % Convert to column vector
        time_gps = time_mjd_gps{i};
        time_grace = time_mjd_grace;
        
        % Perform statistical comparison
        stats = compareTimeSeries(gps_ts, grace_ts, time_gps, time_grace);
        comparison_stats{i} = stats;
        
        if mod(i, 5) == 0 || i == n_valid_stations
            fprintf('    Completed %d/%d stations\n', i, n_valid_stations);
        end
    end
    
    fprintf('  Step 6 completed in %.2f seconds\n\n', toc(step6_start));
    
catch ME
    error('Step 6 failed: %s', ME.message);
end

%% Step 7: Generate comprehensive results
fprintf('Step 7: Generating results and summary statistics...\n');
step7_start = tic;

try
    % Compile overall statistics
    correlations = zeros(n_valid_stations, 1);
    rmse_values = zeros(n_valid_stations, 1);
    nse_values = zeros(n_valid_stations, 1);
    bias_values = zeros(n_valid_stations, 1);
    amplitude_ratios = zeros(n_valid_stations, 1);
    phase_lags = zeros(n_valid_stations, 1);
    n_common_obs = zeros(n_valid_stations, 1);
    
    for i = 1:n_valid_stations
        stats = comparison_stats{i};
        correlations(i) = stats.correlation;
        rmse_values(i) = stats.rmse * 1000;  % Convert to mm
        nse_values(i) = stats.nse;
        bias_values(i) = stats.bias * 1000;  % Convert to mm
        amplitude_ratios(i) = stats.amplitude_ratio;
        phase_lags(i) = stats.phase_lag;
        n_common_obs(i) = stats.n_common;
    end
    
    % Summary statistics
    fprintf('\n=== OVERALL RESULTS SUMMARY ===\n');
    fprintf('Analysis period: MJD %.1f - %.1f\n', min(time_mjd_grace), max(time_mjd_grace));
    fprintf('Number of GPS stations: %d\n', n_valid_stations);
    fprintf('Average common observations per station: %.1f\n', mean(n_common_obs));
    fprintf('\nStatistical Performance:\n');
    fprintf('  Correlation coefficient: %.3f ± %.3f\n', mean(correlations), std(correlations));
    fprintf('  RMSE: %.2f ± %.2f mm\n', mean(rmse_values), std(rmse_values));
    fprintf('  Nash-Sutcliffe Efficiency: %.3f ± %.3f\n', mean(nse_values), std(nse_values));
    fprintf('  Bias (GRACE-GPS): %.2f ± %.2f mm\n', mean(bias_values), std(bias_values));
    fprintf('\nSeasonal Analysis:\n');
    fprintf('  Amplitude ratio (GRACE/GPS): %.2f ± %.2f\n', mean(amplitude_ratios), std(amplitude_ratios));
    fprintf('  Phase lag: %.1f ± %.1f days\n', mean(phase_lags), std(phase_lags));
    
    % Quality assessment
    excellent_corr = sum(correlations > 0.75);
    good_corr = sum(correlations > 0.50 & correlations <= 0.75);
    poor_corr = sum(correlations <= 0.50);
    
    fprintf('\nQuality Assessment:\n');
    fprintf('  Excellent correlation (>0.75): %d stations (%.1f%%)\n', ...
            excellent_corr, 100*excellent_corr/n_valid_stations);
    fprintf('  Good correlation (0.50-0.75): %d stations (%.1f%%)\n', ...
            good_corr, 100*good_corr/n_valid_stations);
    fprintf('  Poor correlation (<0.50): %d stations (%.1f%%)\n', ...
            poor_corr, 100*poor_corr/n_valid_stations);
    
    fprintf('  Step 7 completed in %.2f seconds\n\n', toc(step7_start));
    
catch ME
    error('Step 7 failed: %s', ME.message);
end

%% Step 8: Save results and create visualizations
fprintf('Step 8: Saving results and creating visualizations...\n');
step8_start = tic;

try
    % Save comprehensive results
    results = struct();
    results.config = config;
    results.gps_stations.lat = lat_gps;
    results.gps_stations.lon = lon_gps;
    results.gps_stations.n_stations = n_valid_stations;
    results.time_grace = time_mjd_grace;
    results.comparison_stats = comparison_stats;
    results.summary.correlations = correlations;
    results.summary.rmse_mm = rmse_values;
    results.summary.nse = nse_values;
    results.summary.bias_mm = bias_values;
    results.summary.amplitude_ratios = amplitude_ratios;
    results.summary.phase_lags = phase_lags;
    results.summary.n_common_obs = n_common_obs;
    
    % Save main results
    save(fullfile(config.output_dir, 'grace_gps_analysis_results.mat'), 'results', '-v7.3');
    
    % Create summary statistics table
    station_summary = table(lat_gps, lon_gps, correlations, rmse_values, nse_values, ...
                          bias_values, amplitude_ratios, phase_lags, n_common_obs, ...
                          'VariableNames', {'Latitude', 'Longitude', 'Correlation', ...
                          'RMSE_mm', 'NSE', 'Bias_mm', 'Amplitude_Ratio', 'Phase_Lag_days', 'N_Common'});
    
    writetable(station_summary, fullfile(config.output_dir, 'station_comparison_summary.csv'));
    
    % Generate visualizations if requested
    if config.generate_plots
        fprintf('  Creating visualization plots...\n');
        
        % Plot 1: Correlation map
        figure('Position', [100, 100, 800, 600]);
        scatter(lon_gps, lat_gps, 100, correlations, 'filled');
        colorbar;
        title('GPS-GRACE Correlation by Station');
        xlabel('Longitude (deg)');
        ylabel('Latitude (deg)');
        caxis([0, 1]);
        colormap(jet);
        grid on;
        saveas(gcf, fullfile(config.output_dir, 'correlation_map.png'));
        
        % Plot 2: RMSE distribution
        figure('Position', [100, 200, 800, 400]);
        histogram(rmse_values, 20);
        title('Distribution of RMSE Values');
        xlabel('RMSE (mm)');
        ylabel('Number of Stations');
        grid on;
        saveas(gcf, fullfile(config.output_dir, 'rmse_distribution.png'));
        
        % Plot 3: Seasonal amplitude comparison
        figure('Position', [100, 300, 800, 400]);
        scatter(amplitude_ratios, correlations, 100, 'filled');
        xlabel('Amplitude Ratio (GRACE/GPS)');
        ylabel('Correlation Coefficient');
        title('Seasonal Amplitude vs. Correlation');
        grid on;
        saveas(gcf, fullfile(config.output_dir, 'amplitude_vs_correlation.png'));
        
        close all;  % Clean up figures
    end
    
    fprintf('  Step 8 completed in %.2f seconds\n\n', toc(step8_start));
    
catch ME
    error('Step 8 failed: %s', ME.message);
end

%% Analysis completion
total_time = toc(analysis_start);

fprintf('==========================================\n');
fprintf('GRACE-GPS ANALYSIS COMPLETED SUCCESSFULLY\n');
fprintf('==========================================\n');
fprintf('Total analysis time: %.2f seconds (%.1f minutes)\n', total_time, total_time/60);
fprintf('Results saved to: %s\n', config.output_dir);
fprintf('Key outputs:\n');
fprintf('  - grace_gps_analysis_results.mat (comprehensive results)\n');
fprintf('  - station_comparison_summary.csv (summary table)\n');
if config.generate_plots
    fprintf('  - correlation_map.png (spatial correlation)\n');
    fprintf('  - rmse_distribution.png (error statistics)\n');
    fprintf('  - amplitude_vs_correlation.png (seasonal analysis)\n');
end
fprintf('\nAnalysis completed following Wahr et al. (1998) and Fu & Freymueller (2012)\n');
fprintf('==========================================\n');

end