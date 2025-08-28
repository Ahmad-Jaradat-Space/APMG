function main_grace_analysis_1year()
% main_grace_analysis_1year - 1-Year GRACE-GPS analysis with enhanced visualization
%
% PURPOSE:
%   Comprehensive analysis of elastic crustal deformation using GRACE 
%   satellite data and GPS observations for 2005 (1-year analysis)
%   Following Wahr et al. (1998) and Fu & Freymueller (2012) methodologies.
%
% ENHANCED FEATURES:
%   - Restricted to 2005 data for focused analysis
%   - Enhanced visualization suite (8-10 publication-quality plots)
%   - Comprehensive statistical comparison tables
%   - Improved GPS data loading with proper station mapping
%   - CRITICAL FIX: Proper multi-year reference field for seasonal signal preservation
%     (Uses 2004-2009 static reference instead of single-year mean to preserve seasonal cycle)
%
% WORKFLOW:
%   1. Load Love numbers (PREM model)
%   2. Process 2005 GRACE files with C20/degree-1 corrections
%   3. Load and process GPS time series data (5 stations)
%   4. Convert GRACE to vertical deformation grids
%   5. Extract GRACE at GPS station locations
%   6. Comprehensive statistical comparison
%   7. Generate enhanced visualizations and tables
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
fprintf('GRACE-GPS 1-YEAR ANALYSIS (2005)\n');
fprintf('Enhanced Visualization Suite\n');
fprintf('Based on Wahr et al. (1998) & Fu & Freymueller (2012)\n');
fprintf('==========================================\n\n');

% Record analysis start time
analysis_start = tic;

% Add all necessary paths
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));

%% Configuration parameters
config = struct();

% Analysis period (1 year)
config.analysis_year = 2005;               % Target year for analysis

% GRACE processing parameters
config.nmax = 60;                          % Maximum SH degree for GRACE
config.earth_model = 'PREM';               % Earth model for Love numbers
config.grace_dir = 'data/grace';           
config.c20_file = 'data/aux/C20_RL05.txt'; 
config.deg1_file = 'data/aux/deg1_coef.txt'; 

% GPS processing parameters
config.gps_data_dir = 'data/gps';          
config.gps_coords_file = 'data/gps/GPSLatLong.tenv3'; 
config.detrend_method = 'polynomial';      % Polynomial detrending for GPS
config.detrend_degree = 2;                 % Quadratic detrending

% GPS station information (from coordinate file)
config.station_names = {'P056', 'P140', 'P245', 'P308', 'P345'};

% Spatial analysis parameters (California region)
config.lat_range = [35, 42];               % Expanded to cover all GPS stations
config.lon_range = [-123, -118];           % Expanded longitude range
config.grid_res = 0.5;                     % Higher resolution grid

% Output configuration
config.save_intermediate = true;           
config.generate_plots = true;              
config.output_dir = 'output_1year';

% Create output directory
if ~exist(config.output_dir, 'dir')
    mkdir(config.output_dir);
end

fprintf('Analysis Configuration:\n');
fprintf('  Target Year: %d\n', config.analysis_year);
fprintf('  GPS Stations: %d (%s)\n', length(config.station_names), strjoin(config.station_names, ', '));
fprintf('  Analysis Region: %.1f°-%.1f°N, %.1f°-%.1f°W\n', ...
        config.lat_range, -config.lon_range(2), -config.lon_range(1));
fprintf('  Grid Resolution: %.1f degrees\n', config.grid_res);
fprintf('\n');

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

%% Step 2: Process GRACE time series (2005 only)
fprintf('Step 2: Processing GRACE coefficient files for %d...\n', config.analysis_year);
step2_start = tic;

try
    [cnm_ts, snm_ts, time_mjd_grace] = processGRACEfiles_1year(...
        config.grace_dir, config.c20_file, config.deg1_file, config.analysis_year);
    
    n_months = length(time_mjd_grace);
    fprintf('  Processed %d monthly GRACE solutions for %d\n', n_months, config.analysis_year);
    fprintf('  Time range: MJD %.1f to %.1f\n', min(time_mjd_grace), max(time_mjd_grace));
    fprintf('  Step 2 completed in %.2f seconds\n\n', toc(step2_start));
    
    % Save intermediate results
    if config.save_intermediate
        save(fullfile(config.output_dir, 'grace_coefficients_2005.mat'), ...
             'cnm_ts', 'snm_ts', 'time_mjd_grace', '-v7.3');
    end
    
catch ME
    error('Step 2 failed: %s', ME.message);
end

%% Step 3: Load and process GPS data
fprintf('Step 3: Loading GPS time series data...\n');
step3_start = tic;

try
    % Load GPS coordinates manually due to delimiter issues
    fprintf('  Loading GPS coordinates manually...\n');
    
    % Read GPS coordinate file manually
    fid = fopen(config.gps_coords_file, 'r');
    if fid == -1
        error('Cannot open GPS coordinate file: %s', config.gps_coords_file);
    end
    
    station_map = containers.Map();
    line_num = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        line_num = line_num + 1;
        
        if ischar(line) && ~isempty(line)
            % Skip header line
            if line_num == 1 && contains(line, 'sta_id')
                continue;
            end
            
            % Parse line - assuming whitespace separated
            parts = strsplit(strtrim(line));
            if length(parts) >= 3
                station_id = parts{1};
                lat_val = str2double(parts{2});
                lon_val = str2double(parts{3});
                
                if ~isnan(lat_val) && ~isnan(lon_val)
                    station_map(station_id) = [lat_val, lon_val];
                    fprintf('    Loaded %s: %.3f°N, %.3f°W\n', station_id, lat_val, abs(lon_val));
                end
            end
        end
    end
    fclose(fid);
    
    % Initialize GPS data storage
    n_stations = length(config.station_names);
    gps_data = cell(n_stations, 1);
    time_mjd_gps = cell(n_stations, 1);
    lat_gps = zeros(n_stations, 1);
    lon_gps = zeros(n_stations, 1);
    valid_stations = false(n_stations, 1);
    
    fprintf('  Processing %d GPS stations:\n', n_stations);
    
    for i = 1:n_stations
        station_name = config.station_names{i};
        
        try
            % Get station coordinates
            if isKey(station_map, station_name)
                coords = station_map(station_name);
                lat_gps(i) = coords(1);
                lon_gps(i) = coords(2);
            else
                warning('Station %s not found in coordinate file', station_name);
                continue;
            end
            
            % Construct GPS file path
            gps_file = fullfile(config.gps_data_dir, [station_name '.tenv3']);
            
            if exist(gps_file, 'file')
                % Load GPS time series
                gps_struct = load_tenv3(gps_file);
                
                if ~isempty(gps_struct) && isfield(gps_struct, 't') && length(gps_struct.t) > 100
                    % Extract time and vertical component
                    time_gps_all = gps_struct.t;
                    elevation_all = gps_struct.up;
                    
                    % Filter to 2005 data (approximately)
                    time_filter = (time_gps_all >= mjd_from_year(config.analysis_year)) & ...
                                  (time_gps_all < mjd_from_year(config.analysis_year + 1));
                    
                    if sum(time_filter) > 50  % Minimum observations for 2005
                        time_mjd_gps{i} = time_gps_all(time_filter);
                        elevation = elevation_all(time_filter);
                        
                        % Detrend GPS time series
                        if strcmp(config.detrend_method, 'polynomial')
                            detrend_result = fitPolynomial(time_mjd_gps{i}, elevation, config.detrend_degree, 1.0);
                            gps_data{i} = detrend_result.v;  % Residuals after detrending
                        else
                            gps_data{i} = elevation;  % No detrending
                        end
                        
                        valid_stations(i) = true;
                        fprintf('    %s: %.3f°N, %.3f°W - %d observations\n', ...
                                station_name, lat_gps(i), abs(lon_gps(i)), length(gps_data{i}));
                    else
                        fprintf('    %s: Insufficient 2005 data (%d observations)\n', ...
                                station_name, sum(time_filter));
                    end
                else
                    fprintf('    %s: File empty or insufficient data\n', station_name);
                end
            else
                fprintf('    %s: File not found\n', station_name);
            end
            
        catch ME
            warning('Error processing station %s: %s', station_name, ME.message);
        end
    end
    
    % Filter to valid stations only
    lat_gps = lat_gps(valid_stations);
    lon_gps = lon_gps(valid_stations);
    gps_data = gps_data(valid_stations);
    time_mjd_gps = time_mjd_gps(valid_stations);
    config.station_names = config.station_names(valid_stations);
    n_valid_stations = sum(valid_stations);
    
    fprintf('  Successfully processed %d valid GPS stations\n', n_valid_stations);
    fprintf('  Step 3 completed in %.2f seconds\n\n', toc(step3_start));
    
    if n_valid_stations == 0
        error('No valid GPS stations found for analysis');
    end
    
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
    
    % Load or compute multi-year static reference field (Wahr et al. 1998 methodology)
    fprintf('  Loading multi-year static reference field for seasonal analysis...\n');
    
    % Define reference years for static field (literature recommends 5+ years)
    reference_years = [2004, 2005, 2006, 2007, 2008, 2009]; % 6-year reference period
    reference_file = fullfile(config.output_dir, 'static_reference_field_2004_2009.mat');
    
    if exist(reference_file, 'file')
        % Load existing reference field
        fprintf('  Loading existing reference field: %s\n', reference_file);
        ref_data = load(reference_file);
        cnm_static = ref_data.cnm_static;
        snm_static = ref_data.snm_static;
        fprintf('  Loaded %d-year reference field (%d solutions)\n', ...
                length(ref_data.reference_years), ref_data.valid_count);
    else
        % Compute new reference field
        fprintf('  Computing new multi-year reference field for %d-%d...\n', ...
                min(reference_years), max(reference_years));
        fprintf('  This preserves seasonal signal by using static background (Wahr et al. 1998)\n');
        [cnm_static, snm_static] = computeStaticReferenceField(...
            config.grace_dir, config.c20_file, config.deg1_file, reference_years);
        fprintf('  Reference field computed successfully\n');
    end
    
    % Ensure dimensions match current data
    if size(cnm_static, 1) ~= size(cnm_ts, 1) || size(cnm_static, 2) ~= size(cnm_ts, 2)
        error('Reference field dimensions [%d x %d] do not match current data [%d x %d]', ...
              size(cnm_static, 1), size(cnm_static, 2), size(cnm_ts, 1), size(cnm_ts, 2));
    end
    
    % Process each time step with PROPER mean field removal for seasonal analysis
    fprintf('  Computing vertical deformation for %d time steps using coefficient changes:\n', n_months);
    fprintf('  Method: ΔC_nm = C_nm(month) - C_nm(static_reference) [preserves seasonal cycle]\n');
    
    for t = 1:n_months
        % Extract coefficients for this time step
        cnm_abs = cnm_ts(:, :, t);
        snm_abs = snm_ts(:, :, t);
        
        % Remove STATIC reference field to get coefficient changes (ΔC_nm, ΔS_nm)
        % This preserves seasonal variations relative to multi-year background
        cnm = cnm_abs - cnm_static;
        snm = snm_abs - snm_static;
        
        % Convert to vertical deformation using Wahr et al. (1998) with ΔC_nm, ΔS_nm
        u_vertical = graceToVerticalDeformation(cnm, snm, theta_grid, lambda_grid, h_n, k_n);
        
        u_vertical_ts(:, :, t) = u_vertical;
        
        if mod(t, 6) == 0 || t == n_months
            fprintf('    Processed %d/%d time steps (%.1f%%)\n', t, n_months, 100*t/n_months);
        end
    end
    
    fprintf('  Step 4 completed in %.2f seconds\n\n', toc(step4_start));
    
    % Save intermediate results
    if config.save_intermediate
        save(fullfile(config.output_dir, 'grace_deformation_2005.mat'), ...
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
        
        fprintf('    %s: r=%.3f, RMSE=%.2fmm, NSE=%.3f\n', ...
                config.station_names{i}, stats.correlation, stats.rmse*1000, stats.nse);
    end
    
    fprintf('  Step 6 completed in %.2f seconds\n\n', toc(step6_start));
    
catch ME
    error('Step 6 failed: %s', ME.message);
end

%% Step 7: Generate comprehensive results and tables
fprintf('Step 7: Generating comprehensive results and tables...\n');
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
    correlation_pvals = zeros(n_valid_stations, 1);
    
    for i = 1:n_valid_stations
        stats = comparison_stats{i};
        correlations(i) = stats.correlation;
        rmse_values(i) = stats.rmse * 1000;  % Convert to mm
        nse_values(i) = stats.nse;
        bias_values(i) = stats.bias * 1000;  % Convert to mm
        amplitude_ratios(i) = stats.amplitude_ratio;
        phase_lags(i) = stats.phase_lag;
        n_common_obs(i) = stats.n_common;
        correlation_pvals(i) = stats.correlation_p;
    end
    
    % Create comprehensive comparison table
    station_summary = table(config.station_names', lat_gps, lon_gps, correlations, ...
                          correlation_pvals, rmse_values, nse_values, bias_values, ...
                          amplitude_ratios, phase_lags, n_common_obs, ...
                          'VariableNames', {'Station', 'Latitude', 'Longitude', ...
                          'Correlation', 'P_Value', 'RMSE_mm', 'NSE', 'Bias_mm', ...
                          'Amplitude_Ratio', 'Phase_Lag_days', 'N_Common'});
    
    % Display summary
    fprintf('\n=== COMPREHENSIVE RESULTS SUMMARY (2005) ===\n');
    fprintf('Number of GPS stations: %d\n', n_valid_stations);
    fprintf('GRACE time steps: %d\n', n_months);
    fprintf('Average common observations per station: %.1f\n', mean(n_common_obs));
    
    fprintf('\nSTATISTICAL PERFORMANCE:\n');
    fprintf('  Correlation coefficient: %.3f ± %.3f\n', mean(correlations), std(correlations));
    fprintf('  RMSE: %.2f ± %.2f mm\n', mean(rmse_values), std(rmse_values));
    fprintf('  Nash-Sutcliffe Efficiency: %.3f ± %.3f\n', mean(nse_values), std(nse_values));
    fprintf('  Bias (GRACE-GPS): %.2f ± %.2f mm\n', mean(bias_values), std(bias_values));
    
    fprintf('\nSEASONAL ANALYSIS:\n');
    fprintf('  Amplitude ratio (GRACE/GPS): %.2f ± %.2f\n', mean(amplitude_ratios), std(amplitude_ratios));
    fprintf('  Phase lag: %.1f ± %.1f days\n', mean(phase_lags), std(phase_lags));
    
    % Quality assessment
    excellent_corr = sum(correlations > 0.75);
    good_corr = sum(correlations > 0.50 & correlations <= 0.75);
    poor_corr = sum(correlations <= 0.50);
    
    fprintf('\nQUALITY ASSESSMENT:\n');
    fprintf('  Excellent correlation (>0.75): %d stations (%.1f%%)\n', ...
            excellent_corr, 100*excellent_corr/n_valid_stations);
    fprintf('  Good correlation (0.50-0.75): %d stations (%.1f%%)\n', ...
            good_corr, 100*good_corr/n_valid_stations);
    fprintf('  Poor correlation (<0.50): %d stations (%.1f%%)\n', ...
            poor_corr, 100*poor_corr/n_valid_stations);
    
    % Save comprehensive table
    writetable(station_summary, fullfile(config.output_dir, 'comprehensive_comparison_2005.csv'));
    
    fprintf('  Step 7 completed in %.2f seconds\n\n', toc(step7_start));
    
catch ME
    error('Step 7 failed: %s', ME.message);
end

%% Step 8: Enhanced visualizations
fprintf('Step 8: Creating enhanced visualization suite...\n');
step8_start = tic;

try
    if config.generate_plots
        % Create enhanced visualization suite
        create_enhanced_plots(config, comparison_stats, station_summary, ...
                            u_vertical_ts, lat_grid, lon_grid, time_mjd_grace);
    end
    
    fprintf('  Step 8 completed in %.2f seconds\n\n', toc(step8_start));
    
catch ME
    error('Step 8 failed: %s', ME.message);
end

%% Save comprehensive results
fprintf('Saving comprehensive results...\n');

try
    % Save comprehensive results
    results = struct();
    results.config = config;
    results.analysis_year = config.analysis_year;
    results.gps_stations.names = config.station_names;
    results.gps_stations.lat = lat_gps;
    results.gps_stations.lon = lon_gps;
    results.gps_stations.n_stations = n_valid_stations;
    results.time_grace = time_mjd_grace;
    results.comparison_stats = comparison_stats;
    results.summary = station_summary;
    results.performance.mean_correlation = mean(correlations);
    results.performance.mean_rmse_mm = mean(rmse_values);
    results.performance.mean_nse = mean(nse_values);
    results.performance.mean_bias_mm = mean(bias_values);
    
    % Save main results
    save(fullfile(config.output_dir, 'grace_gps_analysis_2005.mat'), 'results', '-v7.3');
    
catch ME
    error('Results saving failed: %s', ME.message);
end

%% Analysis completion
total_time = toc(analysis_start);

fprintf('==========================================\n');
fprintf('1-YEAR GRACE-GPS ANALYSIS COMPLETED\n');
fprintf('==========================================\n');
fprintf('Analysis Year: %d\n', config.analysis_year);
fprintf('Total analysis time: %.2f seconds (%.1f minutes)\n', total_time, total_time/60);
fprintf('Valid GPS stations: %d\n', n_valid_stations);
fprintf('GRACE time steps: %d\n', n_months);
fprintf('\nResults saved to: %s\n', config.output_dir);
fprintf('Key outputs:\n');
fprintf('  - grace_gps_analysis_2005.mat (comprehensive results)\n');
fprintf('  - comprehensive_comparison_2005.csv (detailed table)\n');
if config.generate_plots
    fprintf('  - Enhanced visualization suite (8+ plots)\n');
end
fprintf('\nMean Performance Metrics:\n');
fprintf('  Correlation: %.3f\n', mean(correlations));
fprintf('  RMSE: %.2f mm\n', mean(rmse_values));
fprintf('  NSE: %.3f\n', mean(nse_values));
fprintf('==========================================\n');

end

%% Helper Functions

function mjd_val = mjd_from_year(year)
% Convert year to MJD (approximate start of year)
    mjd_val = (year - 1858) * 365.25;
end