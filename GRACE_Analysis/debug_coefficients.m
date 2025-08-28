function debug_coefficients()
close all; clear; clc;
addpath(fullfile(pwd, 'functions'));
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));
addpath(fullfile(pwd, 'lib', 'statistics'));
addpath(fullfile(pwd, 'lib'));

fprintf('Step 1: Loading GRACE coefficients...\n');
[cnm_ts, snm_ts, time_mjd_grace] = processGRACEfiles('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt');
n_months = size(cnm_ts, 3);

fprintf('GRACE coefficient statistics:\n');
fprintf('  Matrix size: [%d, %d, %d]\n', size(cnm_ts, 1), size(cnm_ts, 2), size(cnm_ts, 3));
fprintf('  Number of months: %d\n', n_months);
fprintf('  cnm_ts range: [%.2e, %.2e]\n', min(cnm_ts(:)), max(cnm_ts(:)));
fprintf('  snm_ts range: [%.2e, %.2e]\n', min(snm_ts(:)), max(snm_ts(:)));

fprintf('\nSample coefficient values (first time step):\n');
fprintf('  C(2,0): %.6e\n', cnm_ts(3, 1, 1));
fprintf('  C(2,1): %.6e\n', cnm_ts(3, 2, 1));
fprintf('  C(2,2): %.6e\n', cnm_ts(3, 3, 1));
fprintf('  S(2,1): %.6e\n', snm_ts(3, 2, 1));
fprintf('  S(2,2): %.6e\n', snm_ts(3, 3, 1));

fprintf('\nStep 2: Computing static reference field...\n');
reference_years = [2004, 2005, 2006, 2007, 2008, 2009];
[cnm_static, snm_static] = computeStaticReferenceField('data/grace', 'data/aux/C20_RL05.txt', 'data/aux/deg1_coef.txt', reference_years);

fprintf('Static field statistics:\n');
fprintf('  cnm_static range: [%.2e, %.2e]\n', min(cnm_static(:)), max(cnm_static(:)));
fprintf('  snm_static range: [%.2e, %.2e]\n', min(snm_static(:)), max(snm_static(:)));

fprintf('\nSample static field values:\n');
fprintf('  C_static(2,0): %.6e\n', cnm_static(3, 1));
fprintf('  C_static(2,1): %.6e\n', cnm_static(3, 2));
fprintf('  C_static(2,2): %.6e\n', cnm_static(3, 3));

fprintf('\nStep 3: Computing anomaly coefficients...\n');
cnm_anom = cnm_ts(:, :, 1) - cnm_static;
snm_anom = snm_ts(:, :, 1) - snm_static;

fprintf('Anomaly coefficient statistics (first month):\n');
fprintf('  cnm_anom range: [%.2e, %.2e]\n', min(cnm_anom(:)), max(cnm_anom(:)));
fprintf('  snm_anom range: [%.2e, %.2e]\n', min(snm_anom(:)), max(snm_anom(:)));

fprintf('\nSample anomaly values:\n');
fprintf('  C_anom(2,0): %.6e\n', cnm_anom(3, 1));
fprintf('  C_anom(2,1): %.6e\n', cnm_anom(3, 2));
fprintf('  C_anom(2,2): %.6e\n', cnm_anom(3, 3));
fprintf('  S_anom(2,1): %.6e\n', snm_anom(3, 2));
fprintf('  S_anom(2,2): %.6e\n', snm_anom(3, 3));

end