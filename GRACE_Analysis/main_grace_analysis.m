%% Main GRACE and GPS Vertical deformation analysis
close all; clear; clc;
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib'));

%% Confiduration
nmax = 60;
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

%% Phase 1 - Grace Analysis - load love numver, Cnm, Snm, Vertical deformation

fprintf('Step 1: Loading Love numbers (%s model)...\n', earth_model);
[h_n, ~, k_n] = loadLoveNumbers(nmax);

fprintf('Step 2: Processing GRACE coefficient files...\n');
[cnm_ts, snm_ts, time_mjd_grace, grace_start_mjd, grace_end_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file, nmax);
n_months = length(time_mjd_grace);


fprintf('Step 3: Converting GRACE to vertical deformation...\n');

lat_vec = lat_range(1):grid_res:lat_range(2);
lon_vec = lon_range(1):grid_res:lon_range(2);
[lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);

theta_grid = (90 - lat_grid) * pi / 180;
lambda_grid = lon_grid * pi / 180;

[nlat, nlon] = size(lat_grid);
u_vertical_ts = zeros(nlat, nlon, n_months);


for t = 1:n_months
    cnm_abs = cnm_ts(:, :, t);
    snm_abs = snm_ts(:, :, t);
    u_vertical = graceToVerticalDeformation(cnm_abs, snm_abs, theta_grid, lambda_grid, h_n, k_n);
    u_vertical_ts(:, :, t) = u_vertical; % meters

end

%% GPS Anslysis - Load GPS
fprintf('Step 4: Loading GPS time series data...\n');
[lat_gps, lon_gps, gps_data, time_mjd_gps, station_names] = readGPS(gps_coords_file, gps_data_dir);
n_valid_stations = length(lat_gps);


%% Comparision

fprintf('Step 5: Extracting GRACE deformation at GPS stations...\n');

grace_at_gps_ts = extractGRACEatGPS(u_vertical_ts, lat_grid, lon_grid, lat_gps, lon_gps); % meters

% Fit linear trends to GRACE time series for consistent GPS-GRACE comparison
grace_detrended = zeros(n_valid_stations, n_months);

for i = 1:n_valid_stations
    grace_ts = grace_at_gps_ts(i, :)'; % meters (monthly)
    
    % Fit linear trend to GRACE data using time in MJD
    poly_result = fitPolynomial(time_mjd_grace, grace_ts, 1, 0.001);
    
    grace_detrended(i, :) = poly_result.v';      % Store detrended GRACE
end


%% Plotting - Time series

fprintf('Step 6: Creating vertical deformation time series plots...\n');
time_years_grace = mjd2decyear(time_mjd_grace);
figure('Position', [50, 50, 1400, 800]);
set(gcf, 'Color', 'white');

% Station 1
subplot(2, 3, 1);
time_years_gps = mjd2decyear(time_mjd_gps{1});
gps_mm = gps_data{1} * 1000;
plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
hold on;
grace_mm = grace_detrended(1, :) * 1000;
plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{1}), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Year', 'FontSize', 10);
ylabel('Vertical Deformation (mm)', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 9);
grid on; grid minor;
set(gca, 'FontSize', 9);
ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
ylim([-ylim_range, ylim_range]);
hold off;

% Station 2
subplot(2, 3, 2);
time_years_gps = mjd2decyear(time_mjd_gps{2});
gps_mm = gps_data{2} * 1000;
plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
hold on;
grace_mm = grace_detrended(2, :) * 1000;
plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{2}), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Year', 'FontSize', 10);
ylabel('Vertical Deformation (mm)', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 9);
grid on; grid minor;
set(gca, 'FontSize', 9);
ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
ylim([-ylim_range, ylim_range]);
hold off;

% Station 3
subplot(2, 3, 3);
time_years_gps = mjd2decyear(time_mjd_gps{3});
gps_mm = gps_data{3} * 1000;
plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
hold on;
grace_mm = grace_detrended(3, :) * 1000;
plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{3}), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Year', 'FontSize', 10);
ylabel('Vertical Deformation (mm)', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 9);
grid on; grid minor;
set(gca, 'FontSize', 9);
ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
ylim([-ylim_range, ylim_range]);
hold off;

% Station 4
subplot(2, 3, 4);
time_years_gps = mjd2decyear(time_mjd_gps{4});
gps_mm = gps_data{4} * 1000;
plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
hold on;
grace_mm = grace_detrended(4, :) * 1000;
plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{4}), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Year', 'FontSize', 10);
ylabel('Vertical Deformation (mm)', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 9);
grid on; grid minor;
set(gca, 'FontSize', 9);
ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
ylim([-ylim_range, ylim_range]);
hold off;

% Station 5
subplot(2, 3, 5);
time_years_gps = mjd2decyear(time_mjd_gps{5});
gps_mm = gps_data{5} * 1000;
plot(time_years_gps, gps_mm, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
hold on;
grace_mm = grace_detrended(5, :) * 1000;
plot(time_years_grace, grace_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'GRACE (PREM)');
title(sprintf('Station %s (%.2f°N, %.2f°W)', station_names{5}), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Year', 'FontSize', 10);
ylabel('Vertical Deformation (mm)', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 9);
grid on; grid minor;
set(gca, 'FontSize', 9);
ylim_range = max(abs([gps_mm; grace_mm'])) * 1.2;
ylim([-ylim_range, ylim_range]);
hold off;

sgtitle('Vertical Crustal Deformation: GPS vs GRACE (PREM Model)', 'FontSize', 16, 'FontWeight', 'bold');
print(fullfile(output_dir, 'time_series_gps_grace_prem.png'), '-dpng', '-r300');
