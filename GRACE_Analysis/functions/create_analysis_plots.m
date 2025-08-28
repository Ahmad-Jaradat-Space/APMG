function create_analysis_plots(config, station_summary, gps_data, grace_at_gps, ...
                              time_mjd_grace, u_vertical_ts, lat_grid, lon_grid, ...
                              lat_gps, lon_gps)
% create_analysis_plots - Create publication-quality plots for GPS-GRACE comparison
%
% SYNTAX:
%   create_analysis_plots(config, station_summary, gps_data, grace_at_gps, ...)
%
% INPUT:
%   config           - Configuration structure
%   station_summary  - Table with comparison statistics
%   gps_data         - Cell array of GPS time series
%   grace_at_gps     - GRACE values at GPS locations [n_stations x n_time]
%   time_mjd_grace   - Time vector in MJD
%   u_vertical_ts    - GRACE deformation field [nlat x nlon x ntime]
%   lat_grid, lon_grid - Spatial grids
%   lat_gps, lon_gps - GPS coordinates
%
% OUTPUT:
%   Saves publication-quality plots to output directory
%
% Author: GRACE Analysis Project
% Date: 2025

fprintf('  Creating comprehensive visualization suite...\n');

% Set publication-quality defaults
set(0, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesLineWidth', 1.2);
set(0, 'DefaultLineLineWidth', 1.5);

% Extract data
n_stations = height(station_summary);
correlations = station_summary.Correlation;
rmse_values = station_summary.RMSE_mm;
bias_values = station_summary.Bias_mm;
amplitude_ratios = station_summary.Amplitude_Ratio;

% Convert time for plotting
time_years = 1858 + time_mjd_grace / 365.25;

%% Plot 1: Time Series Comparison
fprintf('    Plot 1: Time series comparison...\n');

figure('Position', [100, 100, 1400, 800], 'Name', 'GPS vs GRACE Time Series');

for i = 1:n_stations
    subplot(n_stations, 1, i);
    
    % Plot GPS and GRACE
    plot(time_years, gps_data{i} * 1000, 'b-', 'LineWidth', 2, 'DisplayName', 'GPS');
    hold on;
    plot(time_years, grace_at_gps(i, :) * 1000, 'r--', 'LineWidth', 2, 'DisplayName', 'GRACE');
    
    % Add statistics
    text(0.02, 0.98, sprintf('%s: r=%.3f, RMSE=%.1fmm', ...
                            station_summary.Station{i}, correlations(i), rmse_values(i)), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    ylabel('Displacement (mm)');
    title(sprintf('Station %s (%.3f°N, %.3f°W)', station_summary.Station{i}, ...
                  lat_gps(i), abs(lon_gps(i))));
    legend('Location', 'best');
    grid on;
    
    if i == n_stations
        xlabel('Year');
    end
end

sgtitle(sprintf('GPS vs GRACE Vertical Displacement Comparison (%d)', config.analysis_year), ...
        'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'time_series_comparison.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'time_series_comparison.fig'), 'fig');

%% Plot 2: Correlation Analysis
fprintf('    Plot 2: Correlation and scatter analysis...\n');

figure('Position', [200, 200, 1200, 800], 'Name', 'Correlation Analysis');

% Subplot 1: Scatter plots
subplot(2, 2, 1);
colors = lines(n_stations);
hold on;

for i = 1:n_stations
    gps_mm = gps_data{i} * 1000;
    grace_mm = grace_at_gps(i, :)' * 1000;
    scatter(gps_mm, grace_mm, 50, colors(i,:), 'filled', ...
            'DisplayName', sprintf('%s (r=%.2f)', station_summary.Station{i}, correlations(i)));
end

% Add 1:1 line
lims = [min([ylim, xlim]), max([ylim, xlim])];
plot(lims, lims, 'k--', 'LineWidth', 2, 'DisplayName', '1:1 Line');

xlabel('GPS Displacement (mm)');
ylabel('GRACE Displacement (mm)');
title('GPS vs GRACE Scatter Plot');
legend('Location', 'best');
grid on;
axis equal;

% Subplot 2: Correlation by station
subplot(2, 2, 2);
bar(1:n_stations, correlations, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black');
set(gca, 'XTick', 1:n_stations, 'XTickLabel', station_summary.Station);
ylabel('Correlation Coefficient');
title('Correlation by Station');
ylim([0, 1]);
grid on;

% Add correlation values as text
for i = 1:n_stations
    text(i, correlations(i) + 0.05, sprintf('%.3f', correlations(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Subplot 3: RMSE distribution
subplot(2, 2, 3);
bar(1:n_stations, rmse_values, 'FaceColor', [0.9, 0.6, 0.3], 'EdgeColor', 'black');
set(gca, 'XTick', 1:n_stations, 'XTickLabel', station_summary.Station);
ylabel('RMSE (mm)');
title('RMSE by Station');
grid on;

% Add RMSE values as text
for i = 1:n_stations
    text(i, rmse_values(i) + 0.5, sprintf('%.1f', rmse_values(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Subplot 4: Amplitude comparison
subplot(2, 2, 4);
x_pos = 1:n_stations;
width = 0.35;

bar(x_pos - width/2, station_summary.GPS_Amplitude_mm, width, ...
    'FaceColor', [0.4, 0.8, 0.4], 'EdgeColor', 'black', 'DisplayName', 'GPS');
hold on;
bar(x_pos + width/2, station_summary.GRACE_Amplitude_mm, width, ...
    'FaceColor', [0.8, 0.4, 0.4], 'EdgeColor', 'black', 'DisplayName', 'GRACE');

set(gca, 'XTick', 1:n_stations, 'XTickLabel', station_summary.Station);
ylabel('Amplitude (mm)');
title('Seasonal Amplitude Comparison');
legend('Location', 'best');
grid on;

sgtitle(sprintf('Statistical Analysis Summary (%d)', config.analysis_year), ...
        'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'correlation_analysis.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'correlation_analysis.fig'), 'fig');

%% Plot 3: Spatial Distribution
fprintf('    Plot 3: Spatial distribution...\n');

figure('Position', [300, 300, 1200, 500], 'Name', 'Spatial Distribution');

% Subplot 1: Mean GRACE deformation field
subplot(1, 2, 1);
mean_deformation = mean(u_vertical_ts, 3) * 1000; % Convert to mm
contourf(lon_grid, lat_grid, mean_deformation, 20, 'LineStyle', 'none');
colorbar;
colormap(flipud(parula)); % Use flipped parula colormap
hold on;

% Plot GPS stations
scatter(lon_gps, lat_gps, 150, correlations, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Add station labels
for i = 1:n_stations
    text(lon_gps(i)+0.1, lat_gps(i)+0.1, station_summary.Station{i}, ...
         'FontWeight', 'bold', 'FontSize', 10, 'Color', 'white');
end

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title(sprintf('Mean GRACE Deformation (%d)', config.analysis_year));
grid on;

% Subplot 2: Station performance map
subplot(1, 2, 2);
% Create a simple background
contour(lon_grid, lat_grid, mean_deformation, 10, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
hold on;

% Plot stations colored by performance
scatter(lon_gps, lat_gps, 200, correlations, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);
colorbar;
colormap(gca, 'jet');
caxis([0, 1]);

% Add station labels with performance
for i = 1:n_stations
    text(lon_gps(i)+0.05, lat_gps(i)+0.05, ...
         sprintf('%s\nr=%.2f\nRMSE=%.1fmm', station_summary.Station{i}, ...
                 correlations(i), rmse_values(i)), ...
         'FontSize', 9, 'HorizontalAlignment', 'left', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black');
end

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('Station Performance Map');
grid on;

sgtitle(sprintf('Spatial Analysis (%d)', config.analysis_year), ...
        'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'spatial_distribution.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'spatial_distribution.fig'), 'fig');

%% Plot 4: Performance Dashboard
fprintf('    Plot 4: Performance dashboard...\n');

figure('Position', [400, 400, 1200, 800], 'Name', 'Performance Dashboard');

% Overall performance metrics
subplot(2, 3, 1);
performance_labels = {'Correlation', 'RMSE (mm)', 'Bias (mm)'};
performance_means = [mean(correlations), mean(rmse_values), mean(abs(bias_values))];
performance_stds = [std(correlations), std(rmse_values), std(abs(bias_values))];

bar(performance_means, 'FaceColor', [0.5, 0.7, 0.9], 'EdgeColor', 'black');
hold on;
errorbar(1:3, performance_means, performance_stds, 'k', 'LineWidth', 2, 'LineStyle', 'none');

set(gca, 'XTickLabel', performance_labels);
ylabel('Value');
title('Mean Performance Metrics');
grid on;

% Add values as text
for i = 1:3
    text(i, performance_means(i) + performance_stds(i) + 0.1, ...
         sprintf('%.2f±%.2f', performance_means(i), performance_stds(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Quality distribution
subplot(2, 3, 2);
excellent = sum(correlations > 0.75);
good = sum(correlations > 0.50 & correlations <= 0.75);
poor = sum(correlations <= 0.50);

pie([excellent, good, poor], {'Excellent (>0.75)', 'Good (0.50-0.75)', 'Poor (<0.50)'});
title('Correlation Quality Distribution');

% Performance by station
subplot(2, 3, 3);
station_scores = correlations; % Simple performance score
bar(1:n_stations, station_scores, 'FaceColor', [0.3, 0.8, 0.3], 'EdgeColor', 'black');
set(gca, 'XTick', 1:n_stations, 'XTickLabel', station_summary.Station);
ylabel('Correlation');
title('Performance by Station');
ylim([0, 1]);
grid on;

% Residual analysis
subplot(2, 3, 4);
all_residuals = [];
for i = 1:n_stations
    residuals = (grace_at_gps(i, :)' - gps_data{i}) * 1000; % mm
    all_residuals = [all_residuals; residuals]; %#ok<AGROW>
end

histogram(all_residuals, 15, 'FaceColor', [0.8, 0.6, 0.8], 'EdgeColor', 'black');
xlabel('Residual (mm)');
ylabel('Frequency');
title('Residual Distribution');
grid on;

% Add statistics
text(0.7, 0.9, sprintf('Mean: %.1f mm\nStd: %.1f mm', mean(all_residuals), std(all_residuals)), ...
     'Units', 'normalized', 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Amplitude ratio analysis
subplot(2, 3, 5);
bar(1:n_stations, amplitude_ratios, 'FaceColor', [0.9, 0.7, 0.4], 'EdgeColor', 'black');
hold on;
yline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', 'Perfect Agreement');
set(gca, 'XTick', 1:n_stations, 'XTickLabel', station_summary.Station);
ylabel('Amplitude Ratio (GRACE/GPS)');
title('Amplitude Agreement');
legend('Location', 'best');
grid on;

% Summary statistics text
subplot(2, 3, 6);
axis off;
summary_text = {
    sprintf('ANALYSIS SUMMARY (%d)', config.analysis_year);
    '';
    sprintf('Stations Analyzed: %d', n_stations);
    sprintf('Time Steps: %d', length(time_mjd_grace));
    '';
    'PERFORMANCE METRICS:';
    sprintf('Mean Correlation: %.3f ± %.3f', mean(correlations), std(correlations));
    sprintf('Mean RMSE: %.1f ± %.1f mm', mean(rmse_values), std(rmse_values));
    sprintf('Mean Bias: %.1f ± %.1f mm', mean(bias_values), std(bias_values));
    '';
    'QUALITY DISTRIBUTION:';
    sprintf('Excellent: %d stations (%.0f%%)', excellent, 100*excellent/n_stations);
    sprintf('Good: %d stations (%.0f%%)', good, 100*good/n_stations);
    sprintf('Poor: %d stations (%.0f%%)', poor, 100*poor/n_stations);
};

text(0.1, 0.9, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left');

sgtitle(sprintf('GPS-GRACE Analysis Dashboard (%d)', config.analysis_year), ...
        'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'performance_dashboard.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'performance_dashboard.fig'), 'fig');

fprintf('  Comprehensive visualization completed!\n');
fprintf('  Created 4 publication-quality figures saved to: %s\n', config.output_dir);

% Reset defaults
set(0, 'DefaultFigureRenderer', 'opengl');

end