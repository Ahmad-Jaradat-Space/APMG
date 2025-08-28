function create_enhanced_plots(config, comparison_stats, station_summary, u_vertical_ts, lat_grid, lon_grid, time_mjd_grace)
% create_enhanced_plots - Generate comprehensive visualization suite for GRACE-GPS comparison
%
% SYNTAX:
%   create_enhanced_plots(config, comparison_stats, station_summary, u_vertical_ts, lat_grid, lon_grid, time_mjd_grace)
%
% INPUT:
%   config           - Configuration structure with analysis parameters
%   comparison_stats - Cell array of statistics for each station
%   station_summary  - Summary statistics across all stations
%   u_vertical_ts    - GRACE vertical deformation time series [nlat x nlon x ntime]
%   lat_grid         - Latitude grid [nlat x nlon]
%   lon_grid         - Longitude grid [nlat x nlon]  
%   time_mjd_grace   - GRACE time stamps in MJD [ntime x 1]
%
% OUTPUT:
%   Generates and saves publication-quality plots:
%   - Station locations and correlation map
%   - Time series comparison for each station
%   - GRACE deformation evolution maps
%   - Statistical summary plots
%
% Author: GRACE Analysis Project
% Date: 2025

fprintf('Creating enhanced visualization suite...\n');

% Create output directory if it doesn't exist
if ~exist(config.output_dir, 'dir')
    mkdir(config.output_dir);
end

% Set plot parameters for publication quality
set(0, 'DefaultFigurePosition', [100, 100, 1200, 800]);
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontSize', 12);

%% Plot 1: Station locations and GRACE deformation overview
fprintf('  Creating Plot 1: Station locations and deformation map...\n');

figure(1);
clf;

% Plot GRACE deformation (time-averaged)
u_mean = mean(u_vertical_ts, 3) * 1000; % Convert to mm
pcolor(lon_grid, lat_grid, u_mean);
shading interp;
colorbar;
caxis([-2, 2]); % Limit color scale to reasonable range
colormap('parula'); % Use standard MATLAB colormap
title('Mean GRACE Vertical Deformation (2005) with GPS Stations', 'FontSize', 14);
xlabel('Longitude (°)');
ylabel('Latitude (°)');

% Add GPS station locations (manually defined for the analysis)
hold on;

% Define station coordinates for 2005 analysis
station_coords = struct();
station_coords.P056 = [36.027, -119.063];
station_coords.P140 = [38.829, -120.693];
station_coords.P345 = [40.271, -122.271];

station_names = {'P056', 'P140', 'P345'};

% Plot stations with available data
for i = 1:min(length(station_names), length(comparison_stats))
    if isfield(station_coords, station_names{i}) && ~isempty(comparison_stats{i})
        coords = station_coords.(station_names{i});
        scatter(coords(2), coords(1), 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
        text(coords(2), coords(1), sprintf('  %s\n  r=%.2f', station_names{i}, comparison_stats{i}.correlation), ...
             'FontSize', 10, 'Color', 'white', 'FontWeight', 'bold');
    end
end

grid on;
axis equal;
xlim([min(lon_grid(:))-0.5, max(lon_grid(:))+0.5]);
ylim([min(lat_grid(:))-0.5, max(lat_grid(:))+0.5]);

% Save plot
saveas(gcf, fullfile(config.output_dir, 'grace_stations_overview_2005.png'));
fprintf('    Saved: grace_stations_overview_2005.png\n');

%% Plot 2: Time series comparison for each station
fprintf('  Creating Plot 2: Time series comparisons...\n');

figure(2);
clf;

n_valid_stations = length(comparison_stats);
if n_valid_stations > 0
    for i = 1:min(n_valid_stations, 3) % Show up to 3 stations
        if ~isempty(comparison_stats{i})
            subplot(3, 1, i);
            
            % Convert time to decimal years for plotting
            time_years = 1858.0 + comparison_stats{i}.time_common / 365.25;
            
            % Plot GPS and GRACE time series
            plot(time_years, comparison_stats{i}.gps_interp * 1000, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'GPS');
            hold on;
            plot(time_years, comparison_stats{i}.grace_interp * 1000, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'GRACE');
            
            xlabel('Time (years)');
            ylabel('Vertical Displacement (mm)');
            title(sprintf('%s: r=%.3f, RMSE=%.1fmm', station_names{i}, ...
                         comparison_stats{i}.correlation, comparison_stats{i}.rmse*1000));
            legend('Location', 'best');
            grid on;
            
            % Add correlation info
            text(0.02, 0.98, sprintf('n=%d', comparison_stats{i}.n_common), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top');
        end
    end
end

% Save plot
saveas(gcf, fullfile(config.output_dir, 'timeseries_comparison_2005.png'));
fprintf('    Saved: timeseries_comparison_2005.png\n');

%% Plot 3: GRACE deformation evolution (monthly snapshots)
fprintf('  Creating Plot 3: GRACE deformation evolution...\n');

figure(3);
clf;

% Select 4 representative months (seasonal progression)
months_to_plot = [1, 6, 12, 18]; % Early, spring, summer, late
months_to_plot = months_to_plot(months_to_plot <= size(u_vertical_ts, 3));

for i = 1:length(months_to_plot)
    month_idx = months_to_plot(i);
    
    subplot(2, 2, i);
    u_month = u_vertical_ts(:, :, month_idx) * 1000; % Convert to mm
    
    pcolor(lon_grid, lat_grid, u_month);
    shading interp;
    colorbar;
    caxis([-3, 3]); % Consistent color scale
    colormap('parula');
    
    % Convert MJD to readable date
    date_str = datestr(datenum(1858, 11, 17) + time_mjd_grace(month_idx), 'mmm yyyy');
    title(sprintf('GRACE Deformation - %s', date_str));
    xlabel('Longitude (°)');
    ylabel('Latitude (°)');
    
    % Add GPS stations
    hold on;
    for j = 1:length(station_names)
        if isfield(station_coords, station_names{j})
            coords = station_coords.(station_names{j});
            scatter(coords(2), coords(1), 50, 'k', 'filled');
        end
    end
    
    axis equal;
    xlim([min(lon_grid(:)), max(lon_grid(:))]);
    ylim([min(lat_grid(:)), max(lat_grid(:))]);
end

% Save plot
saveas(gcf, fullfile(config.output_dir, 'grace_evolution_2005.png'));
fprintf('    Saved: grace_evolution_2005.png\n');

%% Plot 4: Statistical summary
fprintf('  Creating Plot 4: Statistical summary...\n');

figure(4);
clf;

if n_valid_stations > 0
    % Extract statistics
    correlations = zeros(n_valid_stations, 1);
    rmse_values = zeros(n_valid_stations, 1);
    amp_ratios = zeros(n_valid_stations, 1);
    
    for i = 1:n_valid_stations
        if ~isempty(comparison_stats{i})
            correlations(i) = comparison_stats{i}.correlation;
            rmse_values(i) = comparison_stats{i}.rmse * 1000; % Convert to mm
            amp_ratios(i) = comparison_stats{i}.amplitude_ratio;
        end
    end
    
    % Create multi-panel statistical plot
    subplot(2, 2, 1);
    bar(1:n_valid_stations, correlations);
    set(gca, 'XTickLabel', station_names(1:n_valid_stations));
    ylabel('Correlation Coefficient');
    title('GPS-GRACE Correlation');
    ylim([-1, 1]);
    grid on;
    
    subplot(2, 2, 2);
    bar(1:n_valid_stations, rmse_values);
    set(gca, 'XTickLabel', station_names(1:n_valid_stations));
    ylabel('RMSE (mm)');
    title('Root Mean Square Error');
    grid on;
    
    subplot(2, 2, 3);
    bar(1:n_valid_stations, amp_ratios);
    set(gca, 'XTickLabel', station_names(1:n_valid_stations));
    ylabel('Amplitude Ratio (GRACE/GPS)');
    title('Seasonal Amplitude Ratio');
    ylim([0, 2]);
    grid on;
    
    subplot(2, 2, 4);
    % Summary statistics histogram
    text(0.1, 0.8, sprintf('Analysis Summary (2005):'), 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.7, sprintf('Stations: %d', n_valid_stations), 'FontSize', 11);
    text(0.1, 0.6, sprintf('Mean correlation: %.3f', mean(correlations)), 'FontSize', 11);
    text(0.1, 0.5, sprintf('Mean RMSE: %.1f mm', mean(rmse_values)), 'FontSize', 11);
    text(0.1, 0.4, sprintf('Mean amp. ratio: %.2f', mean(amp_ratios)), 'FontSize', 11);
    text(0.1, 0.3, sprintf('GRACE time steps: %d', size(u_vertical_ts, 3)), 'FontSize', 11);
    axis off;
end

% Save plot
saveas(gcf, fullfile(config.output_dir, 'statistical_summary_2005.png'));
fprintf('    Saved: statistical_summary_2005.png\n');

fprintf('  Enhanced visualization suite complete!\n');
fprintf('  All plots saved to: %s\n', config.output_dir);

end