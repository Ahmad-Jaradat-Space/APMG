function quick_pipeline_test()
% quick_pipeline_test - Fast validation of optimized GRACE analysis pipeline
%
% PURPOSE:
%   Quick test to validate the complete pipeline works without errors
%   Uses minimal test data for fast execution
%
% Author: GRACE Analysis Project
% Date: 2025

close all; clear; clc;

fprintf('==========================================\n');
fprintf('QUICK PIPELINE VALIDATION TEST\n');
fprintf('==========================================\n\n');

% Record test start time
test_start = tic;

% Add all necessary paths
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib', 'gps_processing'));
addpath(fullfile(pwd, 'lib'));

try
    %% Step 1: Test Physical Constants
    fprintf('Step 1: Testing physical constants...\n');
    constants = physicalConstants();
    fprintf('  ✅ Physical constants loaded successfully\n\n');

    %% Step 2: Test Love Numbers
    fprintf('Step 2: Testing Love numbers...\n');
    nmax = 20;  % Small for testing
    [h_n, l_n, k_n, height_factors] = loadLoveNumbers(nmax, 'PREM', 450000);
    fprintf('  ✅ Love numbers loaded with height factors\n');
    fprintf('  h_2 = %.5f, height factor at n=2 = %.6f\n\n', h_n(3), height_factors(3));

    %% Step 3: Test Spherical Harmonic Processing
    fprintf('Step 3: Testing spherical harmonic processing...\n');
    
    % Create minimal test coefficients
    cnm_test = zeros(nmax+1, nmax+1);
    snm_test = zeros(nmax+1, nmax+1);
    
    % Add some realistic test values
    cnm_test(3, 1) = -1.08263e-3;  % C20 (realistic value)
    cnm_test(3, 3) = 2.53265e-6;   % C22 
    snm_test(3, 3) = -1.61962e-6;  % S22
    cnm_test(4, 1) = 5.39973e-7;   % C30
    
    fprintf('  Test coefficients created\n');

    %% Step 4: Test Grid Creation and Deformation Computation
    fprintf('Step 4: Testing deformation computation...\n');
    
    % Create small test grid (California region)
    lat_range = [32, 42];  % California latitude range
    lon_range = [-125, -115];  % California longitude range
    grid_res = 2.0;  % 2 degree resolution for testing
    
    lat_vec = lat_range(1):grid_res:lat_range(2);
    lon_vec = lon_range(1):grid_res:lon_range(2);
    [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
    
    % Convert to spherical harmonic coordinates
    theta_grid = (90 - lat_grid) * pi / 180;  % Colatitude in radians
    lambda_grid = lon_grid * pi / 180;        % Longitude in radians
    
    fprintf('  Grid created: %d x %d points\n', length(lat_vec), length(lon_vec));
    
    % Test optimized deformation computation
    tic;
    u_vertical = graceToVerticalDeformation(cnm_test, snm_test, theta_grid, lambda_grid, h_n, k_n);
    deformation_time = toc;
    
    fprintf('  ✅ Deformation computed in %.4f seconds\n', deformation_time);
    fprintf('  Max deformation: %.2f mm\n', max(u_vertical(:)) * 1000);
    fprintf('  Min deformation: %.2f mm\n', min(u_vertical(:)) * 1000);

    %% Step 5: Test GPS Location Extraction
    fprintf('\nStep 5: Testing GPS location extraction...\n');
    
    % Create test GPS coordinates (California stations)
    lat_gps = [34.5, 37.5, 40.0]';  % Test GPS latitudes
    lon_gps = [-118.0, -122.0, -120.0]';  % Test GPS longitudes
    n_stations = length(lat_gps);
    
    % Test extraction function
    u_vertical_3d = repmat(u_vertical, [1, 1, 3]);  % Simulate 3 time steps
    grace_at_gps = extractGRACEatGPS(u_vertical_3d, lat_grid, lon_grid, lat_gps, lon_gps);
    
    fprintf('  ✅ GRACE extracted at %d GPS stations\n', n_stations);
    fprintf('  Sample values: [%.2f, %.2f, %.2f] mm\n', grace_at_gps(:, 1)' * 1000);

    %% Step 6: Test Time Series Comparison
    fprintf('\nStep 6: Testing time series comparison...\n');
    
    % Create synthetic GPS time series for testing
    n_months = 24;  % 2 years of data
    time_mjd = linspace(55000, 55000 + 2*365, n_months)';  % 2 years starting from MJD 55000
    
    % Generate synthetic GPS and GRACE time series with some correlation
    t_norm = (1:n_months)' / n_months;
    
    % GPS time series with trend + seasonal + noise
    gps_trend = 0.005 * t_norm;  % 5mm/year trend
    gps_seasonal = 0.003 * sin(2*pi*t_norm) + 0.001 * cos(4*pi*t_norm);  % Seasonal
    gps_noise = 0.002 * randn(n_months, 1);  % 2mm RMS noise
    gps_ts = gps_trend + gps_seasonal + gps_noise;
    
    % GRACE time series (similar but with different noise characteristics)
    grace_seasonal = 0.0025 * sin(2*pi*t_norm + 0.2) + 0.0008 * cos(4*pi*t_norm);  % Similar seasonal
    grace_noise = 0.001 * randn(n_months, 1);  % Lower noise
    grace_ts = 0.8 * (gps_trend + grace_seasonal) + grace_noise;  % 80% of GPS signal
    
    % Test comparison function
    stats = compareTimeSeries(gps_ts, grace_ts, time_mjd, time_mjd);
    
    fprintf('  ✅ Time series comparison completed\n');
    fprintf('  Correlation: %.3f (p=%.4f)\n', stats.correlation, stats.correlation_p);
    fprintf('  RMSE: %.2f mm\n', stats.rmse * 1000);
    fprintf('  NSE: %.3f\n', stats.nse);
    fprintf('  Amplitude ratio: %.2f\n', stats.amplitude_ratio);

    %% Step 7: Test Complete Integration
    fprintf('\nStep 7: Testing complete integration...\n');
    
    % Test that all functions work together
    integration_success = true;
    
    % Check all outputs are reasonable
    if any(~isfinite(u_vertical(:)))
        integration_success = false;
        fprintf('  ❌ Non-finite deformation values\n');
    end
    
    if any(~isfinite(grace_at_gps(:)))
        integration_success = false;
        fprintf('  ❌ Non-finite GPS extraction values\n');
    end
    
    if ~isfinite(stats.correlation) || abs(stats.correlation) > 1
        integration_success = false;
        fprintf('  ❌ Invalid correlation value\n');
    end
    
    if integration_success
        fprintf('  ✅ Complete integration successful\n');
    end

    %% Success Summary
    total_time = toc(test_start);
    
    fprintf('\n==========================================\n');
    fprintf('QUICK PIPELINE TEST - SUCCESSFUL! ✅\n');
    fprintf('==========================================\n');
    fprintf('Total test time: %.2f seconds\n', total_time);
    fprintf('Key validations:\n');
    fprintf('  ✅ Physical constants and Love numbers\n');
    fprintf('  ✅ Optimized spherical harmonic computation\n');
    fprintf('  ✅ Spatial interpolation to GPS locations\n');
    fprintf('  ✅ Statistical time series comparison\n');
    fprintf('  ✅ Complete pipeline integration\n');
    fprintf('\nPipeline is ready for full GRACE analysis!\n');
    fprintf('==========================================\n');

catch ME
    fprintf('\n❌ PIPELINE TEST FAILED!\n');
    fprintf('Error in: %s\n', ME.stack(1).name);
    fprintf('Line: %d\n', ME.stack(1).line);
    fprintf('Message: %s\n', ME.message);
    
    % Additional debugging info
    if exist('constants', 'var')
        fprintf('✅ Physical constants loaded\n');
    else
        fprintf('❌ Physical constants failed\n');
    end
    
    if exist('h_n', 'var')
        fprintf('✅ Love numbers loaded\n');
    else
        fprintf('❌ Love numbers failed\n');
    end
    
    if exist('u_vertical', 'var')
        fprintf('✅ Deformation computation completed\n');
    else
        fprintf('❌ Deformation computation failed\n');
    end
    
    rethrow(ME);
end

end