function main_grace_analysis()
% main_grace_analysis - Complete GRACE-GPS analysis workflow
%
% PURPOSE:
%   Comprehensive analysis of elastic crustal deformation using GRACE 
%   satellite data and GPS observations following Wahr et al. (1998) 
%   and Fu & Freymueller (2012) methodologies.
%
% CRITICAL FIX: Multi-year static reference field for seasonal signal preservation
%   (Uses 2004-2009 static reference instead of absolute coefficients to preserve seasonal cycle)
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

%% Configuration parameters
config = struct();

% GRACE processing parameters
config.nmax = 96;                          % Maximum SH degree for GRACE BB01 files
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

[h_n, l_n, k_n] = loadLoveNumbers(config.nmax, config.earth_model);

fprintf('  Love numbers loaded successfully\n');
fprintf('  h_2 = %.5f, l_2 = %.5f, k_2 = %.5f\n', h_n(3), l_n(3), k_n(3));
fprintf('  Step 1 completed in %.2f seconds\n\n', toc(step1_start));

%% Step 2: Process GRACE time series
fprintf('Step 2: Processing GRACE coefficient files...\n');
step2_start = tic;

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

%% Step 3: Load and process GPS data
fprintf('Step 3: Loading GPS time series data...\n');
step3_start = tic;

% Load GPS coordinates manually due to header line
fprintf('  Loading GPS coordinates from: %s\n', config.gps_coords_file);

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
                end
            end
        end
    end
    fclose(fid);
    
    % Extract station information 
    station_names = keys(station_map);
    n_stations = length(station_names);
    lat_gps = zeros(n_stations, 1);
    lon_gps = zeros(n_stations, 1);
    
    for i = 1:n_stations
        coords = station_map(station_names{i});
        lat_gps(i) = coords(1);
        lon_gps(i) = coords(2);
    end
    
    fprintf('  Successfully loaded %d GPS stations\n', n_stations);
    
    % Initialize GPS data storage
    gps_data = cell(n_stations, 1);
    time_mjd_gps = cell(n_stations, 1);
    
    % Process each GPS station
    valid_stations = false(n_stations, 1);
    
    for i = 1:n_stations
        % Construct filename using actual station name (e.g., P056.tenv3)
        station_file = fullfile(config.gps_data_dir, ...
                              sprintf('%s.tenv3', station_names{i}));
        
        if exist(station_file, 'file')
            % Load time series using existing function
            gps_struct = load_tenv3(station_file);
            
            if ~isempty(gps_struct) && isfield(gps_struct, 't') && isfield(gps_struct, 'up') && length(gps_struct.t) > 50  % Minimum 50 observations
                % Extract vertical component and time from structure
                time_mjd_gps{i} = gps_struct.t;   % MJD time
                elevation = gps_struct.up;        % Vertical component (meters)
                
                % Detrend GPS time series
                if strcmp(config.detrend_method, 'linear')
                    % Add missing sigma0 parameter (assume 1mm standard deviation)
                    poly_result = fitPolynomial(time_mjd_gps{i}, elevation, 1, 0.001);
                    elevation_detrended = poly_result.v;  % Use residuals (detrended data)
                else
                    elevation_detrended = elevation;  % No detrending
                end
                
                gps_data{i} = elevation_detrended;
                valid_stations(i) = true;
            end
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

%% Step 4: Create analysis grid and convert GRACE to deformation
fprintf('Step 4: Converting GRACE to vertical deformation...\n');
step4_start = tic;

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
    end
    
    fprintf('  Step 4 completed in %.2f seconds\n\n', toc(step4_start));
    
% Save intermediate results
if config.save_intermediate
    save(fullfile(config.output_dir, 'grace_deformation.mat'), ...
         'u_vertical_ts', 'lat_grid', 'lon_grid', 'time_mjd_grace', '-v7.3');
end

%% Step 5: Extract GRACE at GPS locations
fprintf('Step 5: Extracting GRACE deformation at GPS stations...\n');
step5_start = tic;

% Extract GRACE values at GPS locations for all time steps
grace_at_gps_ts = extractGRACEatGPS(u_vertical_ts, lat_grid, lon_grid, lat_gps, lon_gps);

fprintf('  Extracted GRACE values for %d stations and %d time steps\n', ...
        n_valid_stations, n_months);
fprintf('  Step 5 completed in %.2f seconds\n\n', toc(step5_start));

%% Step 6: Compare GPS and GRACE time series
fprintf('Step 6: Statistical comparison of GPS and GRACE time series...\n');
step6_start = tic;

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
end

fprintf('  Step 6 completed in %.2f seconds\n\n', toc(step6_start));

%% Step 7: Create vertical deformation time series plots (GPS vs GRACE - PREM)
fprintf('Step 7: Creating vertical deformation time series plots...\n');
step7_start = tic;

% Convert MJD to decimal years for plotting
time_years_grace = mjd2decyear(time_mjd_grace);

% Create production-level time series plot for each GPS station
figure('Position', [50, 50, 1400, 800]);
set(gcf, 'Color', 'white');

n_cols = min(3, n_valid_stations);
n_rows = ceil(n_valid_stations / n_cols);

for i = 1:n_valid_stations
    subplot(n_rows, n_cols, i);
    
    % Convert GPS time to decimal years
    time_years_gps = mjd2decyear(time_mjd_gps{i});
    
    % Plot GPS data in blue
    gps_mm = gps_data{i} * 1000; % Convert to mm
    plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
    hold on;
    
    % Plot GRACE data in red
    grace_mm = grace_at_gps_ts(i, :) * 1000; % Convert to mm
    plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
    
    % Format subplot
    title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{i}, lat_gps(i), -lon_gps(i)), ...
          'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Year', 'FontSize', 10);
    ylabel('Vertical Deformation (mm)', 'FontSize', 10);
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    grid minor;
    set(gca, 'FontSize', 9);
    
    % Set consistent y-axis limits for comparison
    ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
    ylim([-ylim_range, ylim_range]);
    
    hold off;
end

% Add overall title
sgtitle('Vertical Crustal Deformation: GPS vs GRACE (PREM Model)', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Save high-resolution plot
print(fullfile(config.output_dir, 'time_series_gps_grace_prem.png'), '-dpng', '-r300');
print(fullfile(config.output_dir, 'time_series_gps_grace_prem.eps'), '-depsc', '-r300');

fprintf('  Created time series comparison plots\n');
fprintf('  Step 7 completed in %.2f seconds\n\n', toc(step7_start));

%% Step 8: Create spatial difference maps every 5 years
fprintf('Step 8: Creating spatial difference maps (GPS-GRACE) every 5 years...\n');
step8_start = tic;

% Define years for spatial analysis (every 5 years)
analysis_years = [2005, 2010, 2015];
n_years = length(analysis_years);

% Convert years to MJD for finding closest time steps
analysis_mjd = decyear2mjd(analysis_years);

% Find closest GRACE time indices
grace_indices = zeros(n_years, 1);
for i = 1:n_years
    [~, grace_indices(i)] = min(abs(time_mjd_grace - analysis_mjd(i)));
end

% Create figure with subplots
figure('Position', [100, 100, 1200, 400]);
set(gcf, 'Color', 'white');

for i = 1:n_years
    subplot(1, n_years, i);
    
    % Get GRACE deformation for this year
    grace_field = u_vertical_ts(:, :, grace_indices(i)) * 1000; % Convert to mm
    
    % Create interpolated GPS field for comparison
    gps_field = zeros(size(lat_grid));
    
    % Simple interpolation of GPS data to grid
    for j = 1:n_valid_stations
        % Find closest time in GPS data
        gps_time_year = mjd2decyear(time_mjd_gps{j});
        [~, gps_time_idx] = min(abs(gps_time_year - analysis_years(i)));
        
        % Find closest grid point to GPS station
        [~, lat_idx] = min(abs(lat_vec - lat_gps(j)));
        [~, lon_idx] = min(abs(lon_vec - lon_gps(j)));
        
        % Assign GPS value to nearest grid point
        if gps_time_idx <= length(gps_data{j})
            gps_field(lat_idx, lon_idx) = gps_data{j}(gps_time_idx) * 1000; % mm
        end
    end
    
    % Calculate difference (GPS - GRACE) where GPS data exists
    diff_field = nan(size(lat_grid));
    for j = 1:n_valid_stations
        [~, lat_idx] = min(abs(lat_vec - lat_gps(j)));
        [~, lon_idx] = min(abs(lon_vec - lon_gps(j)));
        
        % Find GPS value for this year
        gps_time_year = mjd2decyear(time_mjd_gps{j});
        [~, gps_time_idx] = min(abs(gps_time_year - analysis_years(i)));
        
        if gps_time_idx <= length(gps_data{j})
            grace_val = grace_field(lat_idx, lon_idx);
            gps_val = gps_data{j}(gps_time_idx) * 1000;
            diff_field(lat_idx, lon_idx) = gps_val - grace_val;
        end
    end
    
    % Create surface plot of GRACE field as background
    surf(lon_grid, lat_grid, grace_field, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    
    % Overlay GPS-GRACE differences as colored circles
    for j = 1:n_valid_stations
        % Find GPS value for this year
        gps_time_year = mjd2decyear(time_mjd_gps{j});
        [~, gps_time_idx] = min(abs(gps_time_year - analysis_years(i)));
        
        if gps_time_idx <= length(gps_data{j})
            [~, lat_idx] = min(abs(lat_vec - lat_gps(j)));
            [~, lon_idx] = min(abs(lon_vec - lon_gps(j)));
            
            grace_val = grace_field(lat_idx, lon_idx);
            gps_val = gps_data{j}(gps_time_idx) * 1000;
            diff_val = gps_val - grace_val;
            
            % Color-code the difference
            if abs(diff_val) < 2
                marker_color = 'g'; % Green for good agreement
            elseif abs(diff_val) < 5
                marker_color = 'y'; % Yellow for moderate difference
            else
                marker_color = 'r'; % Red for large difference
            end
            
            plot3(lon_gps(j), lat_gps(j), max(grace_field(:)) + 1, 'o', ...
                  'MarkerSize', 10, 'MarkerFaceColor', marker_color, ...
                  'MarkerEdgeColor', 'k', 'LineWidth', 2);
            
            % Add difference value as text
            text(lon_gps(j), lat_gps(j), max(grace_field(:)) + 2, ...
                 sprintf('%.1f', diff_val), 'FontSize', 8, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
        end
    end
    
    % Format plot
    title(sprintf('Year %d', analysis_years(i)), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Longitude (°)', 'FontSize', 12);
    ylabel('Latitude (°)', 'FontSize', 12);
    zlabel('Deformation (mm)', 'FontSize', 12);
    
    % Set colormap for GRACE background
    colormap(subplot(1, n_years, i), 'RdBu');
    caxis([-max(abs(grace_field(:))), max(abs(grace_field(:)))]);
    
    % Add colorbar
    cb = colorbar;
    cb.Label.String = 'GRACE Deformation (mm)';
    cb.Label.FontSize = 10;
    
    % Set view angle
    view(2); % Top view
    
    % Add grid
    grid on;
    
    hold off;
end

% Add overall title and legend
sgtitle('Spatial Distribution of Vertical Deformation (GPS vs GRACE)', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Add custom legend for GPS-GRACE difference markers
legend_ax = axes('Position', [0.02, 0.02, 0.15, 0.15], 'Visible', 'off');
hold(legend_ax, 'on');
plot(legend_ax, 0, 3, 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
plot(legend_ax, 0, 2, 'oy', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
plot(legend_ax, 0, 1, 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(legend_ax, 0.5, 3, '|GPS-GRACE| < 2mm', 'FontSize', 8);
text(legend_ax, 0.5, 2, '|GPS-GRACE| < 5mm', 'FontSize', 8);
text(legend_ax, 0.5, 1, '|GPS-GRACE| ≥ 5mm', 'FontSize', 8);
set(legend_ax, 'XLim', [0, 3], 'YLim', [0, 4]);

% Save high-resolution plot
print(fullfile(config.output_dir, 'spatial_differences_5year.png'), '-dpng', '-r300');
print(fullfile(config.output_dir, 'spatial_differences_5year.eps'), '-depsc', '-r300');

fprintf('  Created spatial difference maps for years: %s\n', mat2str(analysis_years));
fprintf('  Step 8 completed in %.2f seconds\n\n', toc(step8_start));

%% Analysis completion
total_time = toc(analysis_start);

end