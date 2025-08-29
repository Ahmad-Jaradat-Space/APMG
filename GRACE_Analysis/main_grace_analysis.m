function main_grace_analysis()
close all; clear; clc;
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib'));

nmax = 96;
earth_model = 'PREM';

% Data paths
grace_dir = 'data/grace';
c20_file = 'data/aux/C20_RL05.txt';
deg1_file = 'data/aux/deg1_coef.txt';

gps_data_dir = 'data/gps';
gps_coords_file = 'data/gps/GPSLatLong.tenv3';

% Spatial domain
lat_range = [30, 50];
lon_range = [-130, -110];

grid_res = 1.0;  % degrees

output_dir = 'output';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Step 1: Loading Love numbers (%s model)...\n', earth_model);
step1_start = tic;
[h_n, ~, k_n] = loadLoveNumbers(nmax, earth_model);
fprintf('  Step 1 completed in %.2f seconds\n\n', toc(step1_start));

fprintf('Step 2: Processing GRACE coefficient files...\n');
step2_start = tic;
[cnm_ts, snm_ts, time_mjd_grace, grace_start_mjd, grace_end_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file);
n_months = length(time_mjd_grace);


fprintf('  Step 2 completed in %.2f seconds\n\n', toc(step2_start));

fprintf('Step 3: Loading GPS time series data...\n');
step3_start = tic;

fid = fopen(gps_coords_file, 'r');
if fid == -1
    return;
end
station_map = containers.Map();
line_num = 0;
while ~feof(fid)
    line = fgetl(fid);
    line_num = line_num + 1;
    if ischar(line) && ~isempty(line)
        if line_num == 1 && contains(line, 'sta_id')
            continue;
        end
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

station_names = keys(station_map);
n_stations = length(station_names);

lat_gps = zeros(n_stations, 1);
lon_gps = zeros(n_stations, 1);
for i = 1:n_stations
    coords = station_map(station_names{i});
    lat_gps(i) = coords(1);
    lon_gps(i) = coords(2);
end

gps_data = cell(n_stations, 1);
time_mjd_gps = cell(n_stations, 1);
valid_stations = false(n_stations, 1);

for i = 1:n_stations
    station_file = fullfile(gps_data_dir, sprintf('%s.tenv3', station_names{i}));
    if exist(station_file, 'file')
        gps_struct = load_tenv3(station_file);
        if ~isempty(gps_struct) && isfield(gps_struct, 't') && isfield(gps_struct, 'up') && length(gps_struct.t) > 50
            time_mjd_gps{i} = gps_struct.t;
            % Process GPS data for each station
            elevation = gps_struct.up;
            
            % Special handling for P056 - two-segment linear detrending
            if strcmp(station_names{i}, 'P056')
                % Convert MJD to decimal years for break point
                time_decyear = mjd2decyear(time_mjd_gps{i});
                break_point_decyear = 2014.073; % From break detection analysis
                
                % Use piecewise linear detrending
                poly_result = fitPiecewiseLinear(time_decyear, elevation, break_point_decyear);
                elevation_detrended = poly_result.v;
                
            else
                % Regular polynomial detrending for other stations
                poly_result = fitPolynomial(time_mjd_gps{i}, elevation, 1, 0.001);
                elevation_detrended = poly_result.v;
            end
            
            % Store detrended GPS data
            gps_data{i} = elevation_detrended;
            valid_stations(i) = true;
        end
    end
end

lat_gps = lat_gps(valid_stations);
lon_gps = lon_gps(valid_stations);
gps_data = gps_data(valid_stations);
time_mjd_gps = time_mjd_gps(valid_stations);
n_valid_stations = sum(valid_stations);

fprintf('  Step 3 completed in %.2f seconds\n\n', toc(step3_start));

fprintf('Step 4: Converting GRACE to vertical deformation...\n');
step4_start = tic;

lat_vec = lat_range(1):grid_res:lat_range(2);
lon_vec = lon_range(1):grid_res:lon_range(2);
[lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);

theta_grid = (90 - lat_grid) * pi / 180;
lambda_grid = lon_grid * pi / 180;

[nlat, nlon] = size(lat_grid);
u_vertical_ts = zeros(nlat, nlon, n_months);

% Use pre-drought reference period to preserve California drought signal
reference_years = 2002:2006;
[cnm_static, snm_static] = computeStaticReferenceField(grace_dir, c20_file, deg1_file, reference_years);

fprintf('  Processing %d time steps with static field removal (2002–2006) ...\n', n_months);
for t = 1:n_months
    cnm_abs = cnm_ts(:, :, t);
    snm_abs = snm_ts(:, :, t);

    cnm_delta = cnm_abs - cnm_static;
    snm_delta = snm_abs - snm_static;

    u_vertical = graceToVerticalDeformation(cnm_delta, snm_delta, theta_grid, lambda_grid, h_n, k_n);
    u_vertical_ts(:, :, t) = u_vertical; % meters

end

fprintf('  Step 4 completed in %.2f seconds\n\n', toc(step4_start));

fprintf('Step 5: Extracting GRACE deformation at GPS stations...\n');
step5_start = tic;

grace_at_gps_ts = extractGRACEatGPS(u_vertical_ts, lat_grid, lon_grid, lat_gps, lon_gps); % meters

% Fit linear trends to GRACE time series for consistent GPS-GRACE comparison
grace_detrended = zeros(n_valid_stations, n_months);

for i = 1:n_valid_stations
    grace_ts = grace_at_gps_ts(i, :)'; % meters (monthly)
    
    % Fit linear trend to GRACE data using time in MJD
    poly_result = fitPolynomial(time_mjd_grace, grace_ts, 1, 0.001);
    
    grace_detrended(i, :) = poly_result.v';      % Store detrended GRACE
end

fprintf('  Step 5 completed in %.2f seconds\n\n', toc(step5_start));

fprintf('Step 6: Statistical comparison of GPS and GRACE time series...\n');
step6_start = tic;
comparison_stats = cell(n_valid_stations, 1);
correlations = zeros(n_valid_stations, 1);


for i = 1:n_valid_stations
    gps_ts = gps_data{i};          % meters (daily, detrended)
    grace_ts = grace_detrended(i, :)'; % meters (monthly, detrended)
    time_gps = time_mjd_gps{i};
    time_grace = time_mjd_grace;
    
    % Use improved correlation method with proper temporal alignment
    if length(gps_ts) > 1 && length(grace_ts) > 1
        stats = compareTimeSeriesImproved(gps_ts, time_gps, grace_ts, ...
                                         grace_start_mjd, grace_end_mjd, time_grace);
        correlations(i) = stats.correlation;
        comparison_stats{i} = stats;
    else
        correlations(i) = NaN;
        comparison_stats{i} = struct();
    end
end


fprintf('  Step 6 completed in %.2f seconds\n\n', toc(step6_start));

fprintf('Step 7: Creating vertical deformation time series plots...\n');
step7_start = tic;
time_years_grace = mjd2decyear(time_mjd_grace);
figure('Position', [50, 50, 1400, 800]);
set(gcf, 'Color', 'white');

n_cols = min(3, n_valid_stations);
n_rows = ceil(n_valid_stations / n_cols);

for i = 1:n_valid_stations
    subplot(n_rows, n_cols, i);
    time_years_gps = mjd2decyear(time_mjd_gps{i});
    gps_mm = gps_data{i} * 1000;       % meters -> mm
    plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
    hold on;
    grace_mm = grace_detrended(i, :) * 1000; % meters -> mm  
    plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
    title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{i}, lat_gps(i), -lon_gps(i)), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Year', 'FontSize', 10);
    ylabel('Vertical Deformation (mm)', 'FontSize', 10);
    legend('Location', 'best', 'FontSize', 9);
    grid on; grid minor;
    set(gca, 'FontSize', 9);
    ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
    ylim([-ylim_range, ylim_range]);
    hold off;
end

sgtitle('Vertical Crustal Deformation: GPS vs GRACE (PREM Model)', 'FontSize', 16, 'FontWeight', 'bold');
print(fullfile(output_dir, 'time_series_gps_grace_prem.png'), '-dpng', '-r300');
print(fullfile(output_dir, 'time_series_gps_grace_prem.eps'), '-depsc', '-r300');

fprintf('  Step 7 completed in %.2f seconds\n\n', toc(step7_start));


fprintf('Step 8: Creating spatial difference maps (GPS-GRACE) every 5 years...\n');
step8_start = tic;
analysis_years = [2005, 2010, 2015];
n_years = length(analysis_years);
analysis_mjd = decyear2mjd(analysis_years);

grace_indices = zeros(n_years, 1);
for i = 1:n_years
    [~, grace_indices(i)] = min(abs(time_mjd_grace - analysis_mjd(i)));
end

figure('Position', [100, 100, 1200, 400]);
set(gcf, 'Color', 'white');
for i = 1:n_years
    subplot(1, n_years, i);
    grace_field = u_vertical_ts(:, :, grace_indices(i)) * 1000; % mm
    gps_field = zeros(size(lat_grid));
    for j = 1:n_valid_stations
        gps_time_year = mjd2decyear(time_mjd_gps{j});
        [~, gps_time_idx] = min(abs(gps_time_year - analysis_years(i)));
        [~, lat_idx] = min(abs(lat_vec - lat_gps(j)));
        [~, lon_idx] = min(abs(lon_vec - lon_gps(j)));
        if gps_time_idx <= length(gps_data{j})
            gps_field(lat_idx, lon_idx) = gps_data{j}(gps_time_idx) * 1000; % mm
        end
    end
    diff_field = nan(size(lat_grid));
    for j = 1:n_valid_stations
        [~, lat_idx] = min(abs(lat_vec - lat_gps(j)));
        [~, lon_idx] = min(abs(lon_vec - lon_gps(j)));
        gps_time_year = mjd2decyear(time_mjd_gps{j});
        [~, gps_time_idx] = min(abs(gps_time_year - analysis_years(i)));
        if gps_time_idx <= length(gps_data{j})
            grace_val = grace_field(lat_idx, lon_idx);
            gps_val = gps_data{j}(gps_time_idx) * 1000; % mm
            diff_field(lat_idx, lon_idx) = gps_val - grace_val;
        end
    end

    surf(lon_grid, lat_grid, grace_field, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    for j = 1:n_valid_stations
        gps_time_year = mjd2decyear(time_mjd_gps{j});
        [~, gps_time_idx] = min(abs(gps_time_year - analysis_years(i)));
        if gps_time_idx <= length(gps_data{j})
            [~, lat_idx] = min(abs(lat_vec - lat_gps(j)));
            [~, lon_idx] = min(abs(lon_vec - lon_gps(j)));
            grace_val = grace_field(lat_idx, lon_idx);
            gps_val = gps_data{j}(gps_time_idx) * 1000; % mm
            diff_val = gps_val - grace_val;
            if abs(diff_val) < 2
                marker_color = 'g';
            elseif abs(diff_val) < 5
                marker_color = 'y';
            else
                marker_color = 'r';
            end
            plot3(lon_gps(j), lat_gps(j), max(grace_field(:)) + 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
            text(lon_gps(j), lat_gps(j), max(grace_field(:)) + 2, sprintf('%.1f', diff_val), 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
    title(sprintf('Year %d', analysis_years(i)), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Longitude (°)', 'FontSize', 12);
    ylabel('Latitude (°)', 'FontSize', 12);
    zlabel('Deformation (mm)', 'FontSize', 12);
    colormap(subplot(1, n_years, i), 'jet');
    caxis([-max(abs(grace_field(:))), max(abs(grace_field(:)))]);
    cb = colorbar; cb.Label.String = 'GRACE Deformation (mm)'; cb.Label.FontSize = 10;
    view(2); grid on; hold off;
end

sgtitle('Spatial Distribution of Vertical Deformation (GPS vs GRACE)', 'FontSize', 16, 'FontWeight', 'bold');
legend_ax = axes('Position', [0.02, 0.02, 0.15, 0.15], 'Visible', 'off');
hold(legend_ax, 'on');
plot(legend_ax, 0, 3, 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
plot(legend_ax, 0, 2, 'oy', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
plot(legend_ax, 0, 1, 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
text(legend_ax, 0.5, 3, '|GPS-GRACE| < 2mm', 'FontSize', 8);
text(legend_ax, 0.5, 2, '|GPS-GRACE| < 5mm', 'FontSize', 8);
text(legend_ax, 0.5, 1, '|GPS-GRACE| ≥ 5mm', 'FontSize', 8);
set(legend_ax, 'XLim', [0, 3], 'YLim', [0, 4]);
print(fullfile(output_dir, 'spatial_differences_5year.png'), '-dpng', '-r300');
print(fullfile(output_dir, 'spatial_differences_5year.eps'), '-depsc', '-r300');

fprintf('  Step 8 completed in %.2f seconds\n\n', toc(step8_start));
end