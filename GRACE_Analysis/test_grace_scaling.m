% Test different GRACE scaling factors to match GPS magnitude
addpath(genpath(pwd));

fprintf('=== GRACE SCALING FACTOR TEST ===\n\n');

% Load current results for comparison
fprintf('Loading current analysis results...\n');
[cnm_ts, snm_ts, time_mjd_grace] = processGRACEfiles('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt');
[h_n, l_n, k_n] = loadLoveNumbers(96, 'PREM');

% Use same reference period
reference_years = 2002:2006;
[cnm_static, snm_static] = computeStaticReferenceField('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt', reference_years);

% Test point (California)
lat = 37.0; lon = -120.0;
theta = (90 - lat) * pi / 180;
lambda = lon * pi / 180;

% Test multiple months to get range
test_months = [10, 30, 50, 70, 90, 110, 130, 150]; % Spread across time series
grace_values = zeros(length(test_months), 1);

fprintf('Testing GRACE deformation for %d time points...\n', length(test_months));

for i = 1:length(test_months)
    t = test_months(i);
    if t <= size(cnm_ts, 3)
        cnm_delta = cnm_ts(:,:,t) - cnm_static;
        snm_delta = snm_ts(:,:,t) - snm_static;
        
        % Calculate with CURRENT implementation (includes DDK3 3.9x)
        u_vert = graceToVerticalDeformation(cnm_delta, snm_delta, theta, lambda, h_n, k_n);
        grace_values(i) = u_vert * 1000; % Convert to mm
    end
end

grace_range = max(grace_values) - min(grace_values);
grace_std = std(grace_values);

fprintf('\n=== CURRENT GRACE PERFORMANCE ===\n');
fprintf('GRACE range: %.2f mm (%.2f to %.2f)\n', grace_range, min(grace_values), max(grace_values));
fprintf('GRACE std: %.2f mm\n', grace_std);

% Compare to typical GPS ranges from our analysis
fprintf('\n=== GPS COMPARISON TARGET ===\n');
fprintf('GPS range (P056 drought): ~150 mm\n');
fprintf('GPS seasonal amplitude: ~50 mm\n');
fprintf('Target GRACE seasonal: 10-50 mm\n');

% Test what additional scaling would be needed
fprintf('\n=== SCALING FACTOR ANALYSIS ===\n');

target_ranges = [20, 30, 50]; % Target mm ranges for GRACE
fprintf('Current DDK3 factor: 3.9x\n\n');

for target = target_ranges
    needed_factor = target / grace_range;
    total_factor = 3.9 * needed_factor;
    
    fprintf('For %d mm target range:\n', target);
    fprintf('  Additional scaling needed: %.1fx\n', needed_factor);
    fprintf('  Total factor (vs no DDK3): %.1fx\n', total_factor);
    fprintf('  Resulting range: %.2f mm\n', grace_range * needed_factor);
    fprintf('\n');
end

% Test coefficient magnitudes to see if we're limited by precision
cnm_delta_test = cnm_ts(:,:,50) - cnm_static;
snm_delta_test = snm_ts(:,:,50) - snm_static;

fprintf('=== COEFFICIENT PRECISION CHECK ===\n');
fprintf('Max |Cnm delta|: %.2e\n', max(abs(cnm_delta_test(:))));
fprintf('Max |Snm delta|: %.2e\n', max(abs(snm_delta_test(:))));
fprintf('Machine epsilon: %.2e\n', eps);
fprintf('Precision ratio: %.1f (should be >> 100 for reliable results)\n', max(abs(cnm_delta_test(:)))/eps);

if max(abs(cnm_delta_test(:)))/eps < 100
    fprintf('\nWARNING: Coefficients may be limited by numerical precision!\n');
    fprintf('Consider using higher precision GRACE data or different processing.\n');
end

fprintf('\n=== RECOMMENDATIONS ===\n');
if grace_range < 5
    fprintf('1. GRACE signal too small - try stronger gain correction\n');
    fprintf('2. Alternative: Use different GRACE product (GFZ/JPL vs CSR)\n');
    fprintf('3. Alternative: Use mascon solutions instead of spherical harmonics\n');
elseif grace_range > 50
    fprintf('1. GRACE signal reasonable - check temporal alignment with GPS\n');
    fprintf('2. Verify Love number values and reference frame\n');
else
    fprintf('1. GRACE magnitude acceptable - focus on correlation improvement\n');
    fprintf('2. Check temporal synchronization between GPS and GRACE\n');
end