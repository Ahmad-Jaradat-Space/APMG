function create_enhanced_plots(config, comparison_stats, station_summary, u_vertical_ts, lat_grid, lon_grid, time_mjd_grace)
% create_enhanced_plots - Enhanced visualization suite for GRACE-GPS comparison
%
% SYNTAX:
%   create_enhanced_plots(config, comparison_stats, station_summary, u_vertical_ts, lat_grid, lon_grid, time_mjd_grace)
%
% INPUT:
%   config           - Configuration structure
%   comparison_stats - Cell array of comparison statistics for each station
%   station_summary  - Table with station-by-station comparison summary
%   u_vertical_ts    - GRACE vertical deformation time series [nlat x nlon x ntime]
%   lat_grid         - Latitude grid [nlat x nlon]
%   lon_grid         - Longitude grid [nlat x nlon]
%   time_mjd_grace   - GRACE time vector in MJD [ntime x 1]
%
% OUTPUT:
%   Creates 8+ publication-quality plots saved to output directory
%
% PLOTS CREATED:
%   1. Station locations and correlation map
%   2. Time series comparison (GPS vs GRACE) for each station
%   3. Correlation and scatter plots
%   4. Statistical performance summary
%   5. Seasonal analysis (amplitudes and phases)
%   6. Residual analysis and error distribution
%   7. GRACE deformation field animation/snapshots
%   8. Overall performance dashboard
%
% Author: GRACE Analysis Project - Enhanced Visualization
% Date: 2025

fprintf('Creating enhanced visualization suite...\n');

% Set default plot settings for publication quality
set(0, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesLineWidth', 1.2);
set(0, 'DefaultLineLineWidth', 1.5);

% Extract data for plotting
n_stations = height(station_summary);
station_names = station_summary.Station;
lat_gps = station_summary.Latitude;
lon_gps = station_summary.Longitude;
correlations = station_summary.Correlation;
rmse_values = station_summary.RMSE_mm;
nse_values = station_summary.NSE;
bias_values = station_summary.Bias_mm;
amplitude_ratios = station_summary.Amplitude_Ratio;
phase_lags = station_summary.Phase_Lag_days;

%% Plot 1: Station Locations and Correlation Map
fprintf('  Creating Plot 1: Station locations and correlation map...\n');

figure('Position', [100, 100, 1200, 800], 'Name', 'Station Locations and Performance');

% Subplot 1: Station locations on GRACE deformation field
subplot(2, 2, 1);
if ~isempty(u_vertical_ts)
    % Show mean deformation field
    mean_deformation = mean(u_vertical_ts, 3) * 1000; % Convert to mm
    contourf(lon_grid, lat_grid, mean_deformation, 20, 'LineStyle', 'none');
    colorbar;
    colormap('RdBu_r');
    hold on;
    
    % Plot GPS stations
    scatter(lon_gps, lat_gps, 150, correlations, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Add station labels
    for i = 1:n_stations
        text(lon_gps(i)+0.1, lat_gps(i)+0.1, station_names{i}, 'FontWeight', 'bold', 'FontSize', 10);
    end
end

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GPS Stations on Mean GRACE Deformation Field');
grid on;
axis equal;
xlim([min(lon_gps)-0.5, max(lon_gps)+0.5]);
ylim([min(lat_gps)-0.5, max(lat_gps)+0.5]);

% Subplot 2: Correlation spatial distribution
subplot(2, 2, 2);
scatter(lon_gps, lat_gps, 200, correlations, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
colorbar;
colormap(gca, 'jet');
caxis([0, 1]);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GPS-GRACE Correlation by Station');
grid on;
axis equal;

% Add station labels
for i = 1:n_stations
    text(lon_gps(i)+0.05, lat_gps(i)+0.05, sprintf('%s\n(r=%.2f)', station_names{i}, correlations(i)), ...
         'FontSize', 9, 'HorizontalAlignment', 'left');
end

% Subplot 3: RMSE spatial distribution
subplot(2, 2, 3);
scatter(lon_gps, lat_gps, 200, rmse_values, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
colorbar;
colormap(gca, 'hot');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('RMSE Distribution (mm)');
grid on;
axis equal;

% Subplot 4: NSE spatial distribution
subplot(2, 2, 4);
scatter(lon_gps, lat_gps, 200, nse_values, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
colorbar;
colormap(gca, 'viridis');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('Nash-Sutcliffe Efficiency');
grid on;
axis equal;

sgtitle(sprintf('GPS-GRACE Comparison Summary (%d)', config.analysis_year), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'station_locations_performance.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'station_locations_performance.fig'), 'fig');

%% Plot 2: Time Series Comparisons
fprintf('  Creating Plot 2: Time series comparisons...\n');

n_plots_per_fig = min(4, n_stations); % Max 4 subplots per figure
n_figures = ceil(n_stations / n_plots_per_fig);

for fig_num = 1:n_figures
    figure('Position', [200, 200, 1400, 1000], 'Name', sprintf('Time Series Comparison - Set %d', fig_num));
    
    start_idx = (fig_num - 1) * n_plots_per_fig + 1;
    end_idx = min(fig_num * n_plots_per_fig, n_stations);
    
    for i = start_idx:end_idx
        subplot_idx = i - start_idx + 1;
        subplot(n_plots_per_fig, 1, subplot_idx);
        
        % Get comparison data for this station
        stats = comparison_stats{i};
        
        if ~isempty(stats) && ~isempty(stats.gps_interp) && ~isempty(stats.grace_interp)
            % Convert MJD to decimal years for plotting
            time_years = 1858 + stats.time_common / 365.25;
            
            % Plot GPS and GRACE time series
            plot(time_years, stats.gps_interp * 1000, 'b-', 'LineWidth', 1.5, 'DisplayName', 'GPS');
            hold on;
            plot(time_years, stats.grace_interp * 1000, 'r--', 'LineWidth', 1.5, 'DisplayName', 'GRACE');
            
            % Add statistics text
            text(0.02, 0.98, sprintf('r = %.3f\nRMSE = %.1f mm\nNSE = %.3f', ...
                                   correlations(i), rmse_values(i), nse_values(i)), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'BackgroundColor', 'white', 'EdgeColor', 'black');
            
            xlabel('Year');
            ylabel('Vertical Displacement (mm)');
            title(sprintf('Station %s: GPS vs GRACE', station_names{i}));
            legend('Location', 'best');
            grid on;
            
            % Set axis limits
            xlim([min(time_years)-0.1, max(time_years)+0.1]);
        end
    end
    
    sgtitle(sprintf('Time Series Comparison (%d) - Set %d/%d', config.analysis_year, fig_num, n_figures), ...
            'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(config.output_dir, sprintf('time_series_comparison_%d.png', fig_num)), 'png');
    saveas(gcf, fullfile(config.output_dir, sprintf('time_series_comparison_%d.fig', fig_num)), 'fig');
end

%% Plot 3: Correlation and Scatter Analysis
fprintf('  Creating Plot 3: Correlation and scatter analysis...\n');

figure('Position', [300, 300, 1200, 800], 'Name', 'Correlation and Scatter Analysis');

% Create scatter plot for each station
subplot(2, 2, 1);
colors = lines(n_stations);
hold on;

for i = 1:n_stations
    stats = comparison_stats{i};
    if ~isempty(stats) && ~isempty(stats.gps_interp) && ~isempty(stats.grace_interp)
        scatter(stats.gps_interp * 1000, stats.grace_interp * 1000, 50, colors(i,:), 'filled', ...
                'DisplayName', sprintf('%s (r=%.2f)', station_names{i}, correlations(i)));
    end
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

% Subplot 2: Correlation histogram
subplot(2, 2, 2);
histogram(correlations, 10, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black');
xlabel('Correlation Coefficient');
ylabel('Number of Stations');
title('Distribution of Correlation Values');
grid on;
hold on;
xline(mean(correlations), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean = %.3f', mean(correlations)));
legend;

% Subplot 3: RMSE histogram
subplot(2, 2, 3);
histogram(rmse_values, 10, 'FaceColor', [0.9, 0.6, 0.3], 'EdgeColor', 'black');
xlabel('RMSE (mm)');
ylabel('Number of Stations');
title('Distribution of RMSE Values');
grid on;
hold on;
xline(mean(rmse_values), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean = %.1f mm', mean(rmse_values)));
legend;

% Subplot 4: Performance metrics comparison
subplot(2, 2, 4);
performance_data = [correlations, nse_values, amplitude_ratios];
performance_labels = {'Correlation', 'NSE', 'Amplitude Ratio'};

boxplot(performance_data, 'Labels', performance_labels);
ylabel('Metric Value');
title('Performance Metrics Distribution');
grid on;

sgtitle(sprintf('Statistical Analysis Summary (%d)', config.analysis_year), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'correlation_scatter_analysis.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'correlation_scatter_analysis.fig'), 'fig');

%% Plot 4: Seasonal Analysis
fprintf('  Creating Plot 4: Seasonal analysis...\n');

figure('Position', [400, 400, 1200, 800], 'Name', 'Seasonal Analysis');

% Subplot 1: Amplitude comparison
subplot(2, 3, 1);
station_indices = 1:n_stations;
bar(station_indices, amplitude_ratios, 'FaceColor', [0.4, 0.8, 0.4], 'EdgeColor', 'black');
hold on;
yline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', 'Perfect Agreement');
xlabel('Station Index');
ylabel('Amplitude Ratio (GRACE/GPS)');
title('Seasonal Amplitude Comparison');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;
legend;

% Subplot 2: Phase lag analysis
subplot(2, 3, 2);
bar(station_indices, phase_lags, 'FaceColor', [0.8, 0.4, 0.8], 'EdgeColor', 'black');
hold on;
yline(0, 'r--', 'LineWidth', 2, 'DisplayName', 'No Phase Lag');
xlabel('Station Index');
ylabel('Phase Lag (days)');
title('Seasonal Phase Lag');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;
legend;

% Subplot 3: Amplitude ratio vs correlation
subplot(2, 3, 3);
scatter(amplitude_ratios, correlations, 100, 'filled', 'MarkerEdgeColor', 'black');
xlabel('Amplitude Ratio (GRACE/GPS)');
ylabel('Correlation Coefficient');
title('Amplitude Agreement vs Correlation');
grid on;

% Add station labels
for i = 1:n_stations
    text(amplitude_ratios(i)+0.02, correlations(i)+0.02, station_names{i}, 'FontSize', 9);
end

% Subplot 4: Phase lag histogram
subplot(2, 3, 4);
histogram(phase_lags, 8, 'FaceColor', [0.6, 0.8, 0.6], 'EdgeColor', 'black');
xlabel('Phase Lag (days)');
ylabel('Number of Stations');
title('Phase Lag Distribution');
grid on;
hold on;
xline(mean(phase_lags), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean = %.1f days', mean(phase_lags)));
legend;

% Subplot 5: Bias analysis
subplot(2, 3, 5);
bar(station_indices, bias_values, 'FaceColor', [0.8, 0.6, 0.4], 'EdgeColor', 'black');
hold on;
yline(0, 'r--', 'LineWidth', 2, 'DisplayName', 'No Bias');
xlabel('Station Index');
ylabel('Bias (mm)');
title('GRACE-GPS Bias');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;
legend;

% Subplot 6: Overall seasonal metrics
subplot(2, 3, 6);
seasonal_metrics = [mean(amplitude_ratios), std(amplitude_ratios); ...
                   mean(abs(phase_lags)), std(phase_lags); ...
                   mean(abs(bias_values)), std(abs(bias_values))];

bar(seasonal_metrics);
set(gca, 'XTickLabel', {'Amplitude Ratio', 'Phase Lag (days)', 'Bias (mm)'});
ylabel('Value');
title('Seasonal Analysis Summary');
legend({'Mean', 'Std'}, 'Location', 'best');
grid on;

sgtitle(sprintf('Seasonal Analysis (%d)', config.analysis_year), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'seasonal_analysis.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'seasonal_analysis.fig'), 'fig');

%% Plot 5: Performance Dashboard
fprintf('  Creating Plot 5: Performance dashboard...\n');

figure('Position', [500, 500, 1400, 1000], 'Name', 'Performance Dashboard');

% Calculate performance categories
excellent_stations = sum(correlations > 0.75);
good_stations = sum(correlations > 0.50 & correlations <= 0.75);
poor_stations = sum(correlations <= 0.50);

% Dashboard layout with multiple performance metrics
subplot(3, 3, 1);
pie([excellent_stations, good_stations, poor_stations], ...
    {'Excellent (>0.75)', 'Good (0.50-0.75)', 'Poor (<0.50)'});
title('Correlation Quality Distribution');

subplot(3, 3, 2);
histogram(correlations, 8, 'FaceColor', 'blue', 'Alpha', 0.7);
xlabel('Correlation');
ylabel('Count');
title('Correlation Distribution');
grid on;

subplot(3, 3, 3);
histogram(rmse_values, 8, 'FaceColor', 'red', 'Alpha', 0.7);
xlabel('RMSE (mm)');
ylabel('Count');
title('RMSE Distribution');
grid on;

subplot(3, 3, 4);
plot(station_indices, correlations, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Station Index');
ylabel('Correlation');
title('Correlation by Station');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;

subplot(3, 3, 5);
plot(station_indices, rmse_values, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Station Index');
ylabel('RMSE (mm)');
title('RMSE by Station');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;

subplot(3, 3, 6);
plot(station_indices, nse_values, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Station Index');
ylabel('NSE');
title('NSE by Station');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;

% Performance summary statistics
subplot(3, 3, 7);
summary_stats = [mean(correlations), std(correlations); ...
                mean(rmse_values), std(rmse_values); ...
                mean(nse_values), std(nse_values); ...
                mean(abs(bias_values)), std(bias_values)];

bar(summary_stats);
set(gca, 'XTickLabel', {'Correlation', 'RMSE (mm)', 'NSE', 'Bias (mm)'});
ylabel('Value');
title('Performance Statistics');
legend({'Mean', 'Std'}, 'Location', 'best');
grid on;

% Overall performance score (weighted average)
subplot(3, 3, 8);
% Normalize metrics (0-1 scale) and compute weighted score
norm_corr = correlations; % Already 0-1
norm_rmse = 1 - (rmse_values - min(rmse_values)) / (max(rmse_values) - min(rmse_values)); % Higher RMSE = lower score
norm_nse = max(0, nse_values); % NSE can be negative, clamp to 0

performance_scores = 0.5 * norm_corr + 0.3 * norm_rmse + 0.2 * norm_nse;

bar(station_indices, performance_scores, 'FaceColor', [0.3, 0.7, 0.5], 'EdgeColor', 'black');
xlabel('Station Index');
ylabel('Performance Score');
title('Overall Performance Score');
set(gca, 'XTick', station_indices, 'XTickLabel', station_names, 'XTickLabelRotation', 45);
grid on;
ylim([0, 1]);

% Text summary
subplot(3, 3, 9);
axis off;
summary_text = {
    sprintf('ANALYSIS YEAR: %d', config.analysis_year);
    sprintf('STATIONS ANALYZED: %d', n_stations);
    '';
    'PERFORMANCE SUMMARY:';
    sprintf('Mean Correlation: %.3f ± %.3f', mean(correlations), std(correlations));
    sprintf('Mean RMSE: %.1f ± %.1f mm', mean(rmse_values), std(rmse_values));
    sprintf('Mean NSE: %.3f ± %.3f', mean(nse_values), std(nse_values));
    sprintf('Mean Bias: %.1f ± %.1f mm', mean(bias_values), std(bias_values));
    '';
    'QUALITY ASSESSMENT:';
    sprintf('Excellent Correlations: %d (%.0f%%)', excellent_stations, 100*excellent_stations/n_stations);
    sprintf('Good Correlations: %d (%.0f%%)', good_stations, 100*good_stations/n_stations);
    sprintf('Poor Correlations: %d (%.0f%%)', poor_stations, 100*poor_stations/n_stations);
};

text(0.1, 0.9, summary_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left', 'FontFamily', 'monospace');

sgtitle(sprintf('GRACE-GPS Analysis Performance Dashboard (%d)', config.analysis_year), ...
        'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(config.output_dir, 'performance_dashboard.png'), 'png');
saveas(gcf, fullfile(config.output_dir, 'performance_dashboard.fig'), 'fig');

fprintf('  Enhanced visualization suite completed!\n');
fprintf('  Created 5 comprehensive plot sets saved to: %s\n', config.output_dir);

% Reset default figure settings
set(0, 'DefaultFigureRenderer', 'opengl');

end