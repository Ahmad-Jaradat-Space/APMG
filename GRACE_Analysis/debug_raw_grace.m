function debug_raw_grace()
close all; clear; clc;
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib'));

fprintf('Loading GRACE data and Love numbers...\n');
[cnm_ts, snm_ts, time_mjd_grace] = processGRACEfiles('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt');
[h_n, l_n, k_n] = loadLoveNumbers(96, 'PREM');

% Set up a simple grid point for testing
lat = 40.0; % degrees
lon = -120.0; % degrees
theta = (90 - lat) * pi / 180;
lambda = lon * pi / 180;

fprintf('\nTesting deformation calculation...\n');
fprintf('Test location: %.1f°N, %.1f°W\n', lat, -lon);

% Test with raw coefficients (no static field subtraction)
fprintf('\n=== RAW GRACE COEFFICIENTS (no static field subtraction) ===\n');
u_raw_ts = zeros(1, size(cnm_ts, 3));
for t = 1:min(12, size(cnm_ts, 3)) % Test first 12 months
    cnm = cnm_ts(:, :, t);
    snm = snm_ts(:, :, t);
    u_raw = graceToVerticalDeformation(cnm, snm, theta, lambda, h_n, k_n);
    u_raw_ts(t) = u_raw;
    fprintf('Month %2d: %.3f mm\n', t, u_raw * 1000);
end

fprintf('\nRaw deformation statistics (first 12 months):\n');
fprintf('  Range: [%.2f, %.2f] mm\n', min(u_raw_ts(1:12))*1000, max(u_raw_ts(1:12))*1000);
fprintf('  Std: %.2f mm\n', std(u_raw_ts(1:12))*1000);

% Test with static field subtracted (current approach)
fprintf('\n=== WITH STATIC FIELD SUBTRACTION (current approach) ===\n');
reference_years = [2004, 2005, 2006, 2007, 2008, 2009];
[cnm_static, snm_static] = computeStaticReferenceField('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt', reference_years);

u_anom_ts = zeros(1, size(cnm_ts, 3));
for t = 1:min(12, size(cnm_ts, 3))
    cnm_anom = cnm_ts(:, :, t) - cnm_static;
    snm_anom = snm_ts(:, :, t) - snm_static;
    u_anom = graceToVerticalDeformation(cnm_anom, snm_anom, theta, lambda, h_n, k_n);
    u_anom_ts(t) = u_anom;
    fprintf('Month %2d: %.3f mm\n', t, u_anom * 1000);
end

fprintf('\nAnomaly deformation statistics (first 12 months):\n');
fprintf('  Range: [%.2f, %.2f] mm\n', min(u_anom_ts(1:12))*1000, max(u_anom_ts(1:12))*1000);
fprintf('  Std: %.2f mm\n', std(u_anom_ts(1:12))*1000);

% Check coefficient temporal variations
fprintf('\n=== COEFFICIENT TEMPORAL VARIATIONS ===\n');
c20_ts = squeeze(cnm_ts(3, 1, 1:12));
c21_ts = squeeze(cnm_ts(3, 2, 1:12));
s21_ts = squeeze(snm_ts(3, 2, 1:12));

fprintf('C(2,0) temporal variation:\n');
fprintf('  Range: [%.2e, %.2e]\n', min(c20_ts), max(c20_ts));
fprintf('  Std: %.2e\n', std(c20_ts));

fprintf('C(2,1) temporal variation:\n');
fprintf('  Range: [%.2e, %.2e]\n', min(c21_ts), max(c21_ts));
fprintf('  Std: %.2e\n', std(c21_ts));

end