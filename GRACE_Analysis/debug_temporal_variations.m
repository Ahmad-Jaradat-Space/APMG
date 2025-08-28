function debug_temporal_variations()
close all; clear; clc;
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib'));

fprintf('Loading GRACE data...\n');
[cnm_ts, snm_ts, time_mjd_grace] = processGRACEfiles('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt');
time_years_grace = mjd2decyear(time_mjd_grace);

fprintf('Analyzing temporal variations in key coefficients...\n');

% Focus on low-degree coefficients (2,0), (2,1), (2,2) which dominate large-scale signals
c20_ts = squeeze(cnm_ts(3, 1, :)); % C(2,0)
c21_ts = squeeze(cnm_ts(3, 2, :)); % C(2,1) 
c22_ts = squeeze(cnm_ts(3, 3, :)); % C(2,2)
s21_ts = squeeze(snm_ts(3, 2, :)); % S(2,1)
s22_ts = squeeze(snm_ts(3, 3, :)); % S(2,2)

fprintf('\nCoefficient statistics:\n');
fprintf('C(2,0): mean=%.6e, std=%.6e, range=[%.6e, %.6e]\n', ...
    mean(c20_ts), std(c20_ts), min(c20_ts), max(c20_ts));
fprintf('C(2,1): mean=%.6e, std=%.6e, range=[%.6e, %.6e]\n', ...
    mean(c21_ts), std(c21_ts), min(c21_ts), max(c21_ts));
fprintf('C(2,2): mean=%.6e, std=%.6e, range=[%.6e, %.6e]\n', ...
    mean(c22_ts), std(c22_ts), min(c22_ts), max(c22_ts));
fprintf('S(2,1): mean=%.6e, std=%.6e, range=[%.6e, %.6e]\n', ...
    mean(s21_ts), std(s21_ts), min(s21_ts), max(s21_ts));
fprintf('S(2,2): mean=%.6e, std=%.6e, range=[%.6e, %.6e]\n', ...
    mean(s22_ts), std(s22_ts), min(s22_ts), max(s22_ts));

% Test deformation calculation with raw coefficients (first few degrees only)
fprintf('\nTesting deformation with different degree truncations:\n');
[h_n, l_n, k_n] = loadLoveNumbers(96, 'PREM');

% Test location: 40°N, 120°W
lat = 40.0;
lon = -120.0;
theta = (90 - lat) * pi / 180;
lambda = lon * pi / 180;

% Test with degree 2 only
nmax_test = 2;
cnm_test = cnm_ts(1:nmax_test+1, 1:nmax_test+1, :);
snm_test = snm_ts(1:nmax_test+1, 1:nmax_test+1, :);
h_test = h_n(1:nmax_test+1);
k_test = k_n(1:nmax_test+1);

u_deg2_ts = zeros(size(cnm_test, 3), 1);
for t = 1:min(12, size(cnm_test, 3))
    cnm = cnm_test(:, :, t);
    snm = snm_test(:, :, t);
    u_deg2 = graceToVerticalDeformation(cnm, snm, theta, lambda, h_test, k_test);
    u_deg2_ts(t) = u_deg2;
end

fprintf('Degree 2 only deformation (first 12 months):\n');
fprintf('  Range: [%.3f, %.3f] mm\n', min(u_deg2_ts(1:12))*1000, max(u_deg2_ts(1:12))*1000);
fprintf('  Std: %.3f mm\n', std(u_deg2_ts(1:12))*1000);

% Compare with temporal mean subtracted
cnm_mean_test = mean(cnm_test, 3);
snm_mean_test = mean(snm_test, 3);

u_anom_deg2_ts = zeros(size(cnm_test, 3), 1);
for t = 1:min(12, size(cnm_test, 3))
    cnm_anom = cnm_test(:, :, t) - cnm_mean_test;
    snm_anom = snm_test(:, :, t) - snm_mean_test;
    u_anom_deg2 = graceToVerticalDeformation(cnm_anom, snm_anom, theta, lambda, h_test, k_test);
    u_anom_deg2_ts(t) = u_anom_deg2;
end

fprintf('\nDegree 2 anomaly deformation (first 12 months):\n');
fprintf('  Range: [%.3f, %.3f] mm\n', min(u_anom_deg2_ts(1:12))*1000, max(u_anom_deg2_ts(1:12))*1000);
fprintf('  Std: %.3f mm\n', std(u_anom_deg2_ts(1:12))*1000);

end