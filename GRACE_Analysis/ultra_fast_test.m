function ultra_fast_test()
% ultra_fast_test - Ultra-fast GRACE pipeline validation (<1 second)
%
% PURPOSE:
%   Minimal, high-speed test to validate core GRACE pipeline functionality
%   Optimized for maximum speed with essential validation only
%
% TARGET: <1 second total execution time
%
% Author: GRACE Analysis Project
% Date: 2025

% Start timing
tic;

% Minimal output - just essential status
fprintf('🚀 ULTRA-FAST GRACE PIPELINE TEST\n');

% Add paths (silent)
addpath(fullfile(pwd, 'functions'), fullfile(pwd, 'lib', 'spherical_harmonics'), ...
        fullfile(pwd, 'lib', 'time_utils'), fullfile(pwd, 'lib', 'statistics'), fullfile(pwd, 'lib'));

try
    %% Test 1: Physical Constants (Target: 0.1s)
    constants = physicalConstants();
    assert(constants.GM > 3.9e14 && constants.R == 6378137, 'Constants validation failed');
    fprintf('✅ Constants ');
    
    %% Test 2: Love Numbers (Target: 0.1s) 
    nmax = 5;  % Minimal degree for speed
    [h_n, l_n, k_n, height_factors] = loadLoveNumbers(nmax, 'PREM', 450000);
    assert(length(h_n) == nmax+1 && abs(h_n(3) + 0.31) < 0.1, 'Love numbers failed');
    fprintf('✅ Love ');
    
    %% Test 3: Minimal Deformation Computation (Target: 0.2s)
    % Tiny 3x3 grid for maximum speed
    lat_vec = [35, 37.5, 40];  % 3 points only
    lon_vec = [-120, -117.5, -115];  % 3 points only
    [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec);
    
    % Convert to spherical coordinates
    theta_grid = (90 - lat_grid) * pi / 180;
    lambda_grid = lon_grid * pi / 180;
    
    % Realistic small test coefficients (similar to actual GRACE variations)
    cnm = zeros(nmax+1, nmax+1);
    snm = zeros(nmax+1, nmax+1);
    cnm(3, 1) = -1.0e-8;  % Small C20 variation (~0.01mm scale)
    
    % Fast deformation computation
    u_vertical = graceToVerticalDeformation(cnm, snm, theta_grid, lambda_grid, h_n, k_n);
    assert(all(isfinite(u_vertical(:))), 'Deformation computation failed');
    fprintf('✅ Deformation ');
    
    %% Test 4: Single GPS Extraction (Target: 0.1s)
    lat_gps = 37.5;  % Single point
    lon_gps = -117.5;
    
    % Single extraction (no time series)
    grace_value = extractGRACEatGPS(u_vertical, lat_grid, lon_grid, lat_gps, lon_gps);
    assert(isfinite(grace_value), 'GPS extraction failed');
    fprintf('✅ GPS ');
    
    %% Test 5: Basic Validation (Target: 0.1s)
    % Check reasonable magnitudes (should be sub-mm scale with small coefficients)
    max_deform_mm = max(abs(u_vertical(:))) * 1000;
    assert(max_deform_mm >= 0.001 && max_deform_mm < 100, 'Unrealistic deformation values');
    assert(abs(grace_value) < 0.01, 'GPS value out of range');  % Within 1cm
    fprintf('✅ Validation\n');
    
    %% Success
    total_time = toc;
    fprintf('🎯 SUCCESS: All tests passed in %.3f seconds\n', total_time);
    
    if total_time < 1.0
        fprintf('⚡ PERFORMANCE: Target <1s ACHIEVED!\n');
    else
        fprintf('⏱️  PERFORMANCE: %.3fs (target was <1s)\n', total_time);
    end
    
    % Minimal stats
    fprintf('📊 Results: Max deformation %.1f mm, GPS value %.2f mm\n', ...
            max_deform_mm, grace_value*1000);
    
catch ME
    fprintf('❌ FAILED: %s\n', ME.message);
    fprintf('   Location: %s:%d\n', ME.stack(1).name, ME.stack(1).line);
    fprintf('   Time: %.3f seconds\n', toc);
    rethrow(ME);
end

fprintf('🏁 Ultra-fast pipeline validation complete!\n');

end