function [cnm_ts, snm_ts, time_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file)
% processGRACEfiles - Process all GRACE files with corrections
%
% SYNTAX:
%   [cnm_ts, snm_ts, time_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file)
%
% INPUT:
%   grace_dir - Directory containing .gfc files
%   c20_file  - Path to C20 replacement file (SLR-based)
%   deg1_file - Path to degree-1 coefficients file
%
% OUTPUT:
%   cnm_ts    - Time series of Cnm coefficients [n+1 x m+1 x n_months]
%   snm_ts    - Time series of Snm coefficients [n+1 x m+1 x n_months]
%   time_mjd  - Time stamps in Modified Julian Date [n_months x 1]
%
% REFERENCES:
%   - Swenson et al. (2008) for degree-1 coefficients
%   - Cheng et al. (2013) for C20 replacement
%   - Wahr et al. (1998) for GRACE processing methodology
%
% NOTES:
%   - Uses existing readSHC.m function
%   - Applies DDK3 filtering (already in files)
%   - Handles BA01/BB01 (A/B satellite) files
%
% Author: GRACE Analysis Project
% Date: 2025

% Add path to spherical harmonics functions
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));

% Input validation
if ~exist(grace_dir, 'dir')
    error('GRACE directory does not exist: %s', grace_dir);
end

if ~exist(c20_file, 'file')
    error('C20 file does not exist: %s', c20_file);
end

if ~exist(deg1_file, 'file')
    error('Degree-1 file does not exist: %s', deg1_file);
end

fprintf('Processing GRACE files from: %s\n', grace_dir);
fprintf('Using C20 corrections from: %s\n', c20_file);
fprintf('Using degree-1 coefficients from: %s\n', deg1_file);

%% Step 1: Get list of GRACE files and extract time information
gfc_files = dir(fullfile(grace_dir, '*.gfc'));
if isempty(gfc_files)
    error('No .gfc files found in directory: %s', grace_dir);
end

fprintf('Found %d GRACE files\n', length(gfc_files));

% Extract time information from filenames
% Format: kfilter_DDK3_GSM-2_YYYYDDD-YYYYDDD_GRAC_UTCSR_BA01_0600.gfc
n_files = length(gfc_files);
file_dates = zeros(n_files, 2); % [start_date, end_date] as decimal years

for i = 1:n_files
    filename = gfc_files(i).name;
    
    % Extract date range from filename using regex
    pattern = 'GSM-2_(\d{7})-(\d{7})_';
    tokens = regexp(filename, pattern, 'tokens');
    
    if ~isempty(tokens)
        start_date = str2double(tokens{1}{1});
        end_date = str2double(tokens{1}{2});
        
        % Convert YYYYDDD to decimal year
        start_year = floor(start_date / 1000);
        start_doy = mod(start_date, 1000);
        end_year = floor(end_date / 1000);
        end_doy = mod(end_date, 1000);
        
        file_dates(i, 1) = start_year + (start_doy - 1) / 365.25;
        file_dates(i, 2) = end_year + (end_doy - 1) / 365.25;
    else
        warning('Could not parse date from filename: %s', filename);
        file_dates(i, :) = NaN;
    end
end

% Remove files with unparseable dates
valid_idx = ~any(isnan(file_dates), 2);
gfc_files = gfc_files(valid_idx);
file_dates = file_dates(valid_idx, :);
n_files = length(gfc_files);

% Calculate mid-point times and convert to MJD
mid_times = mean(file_dates, 2);
time_mjd = zeros(n_files, 1);
for i = 1:n_files
    % Convert decimal year to MJD using existing function
    % Note: decyear2mjd function should be available
    try
        time_mjd(i) = decyear2mjd(mid_times(i));
    catch
        % Fallback calculation if function not available
        time_mjd(i) = (mid_times(i) - 1858.0) * 365.25; % Approximate MJD
    end
end

% Sort files by time
[time_mjd, sort_idx] = sort(time_mjd);
gfc_files = gfc_files(sort_idx);
file_dates = file_dates(sort_idx, :);
mid_times = mid_times(sort_idx);

fprintf('Time range: %.3f - %.3f (decimal years)\n', min(mid_times), max(mid_times));

%% Step 2: Load C20 replacement data
fprintf('Loading C20 replacement data...\n');
c20_data = load(c20_file);
c20_time = c20_data(:, 1); % Decimal year
c20_values = c20_data(:, 2); % Normalized C20 values

%% Step 3: Load degree-1 coefficients
fprintf('Loading degree-1 coefficients...\n');
deg1_data = load(deg1_file);
deg1_time = deg1_data(:, 1); % Decimal year
deg1_c10 = deg1_data(:, 2); % C10 (unnormalized)
deg1_c11 = deg1_data(:, 3); % C11 (unnormalized)
deg1_s11 = deg1_data(:, 4); % S11 (unnormalized)

%% Step 4: Read first file to get dimensions
fprintf('Reading first file to determine array dimensions...\n');
first_file = fullfile(grace_dir, gfc_files(1).name);

% Read using existing readSHC function
try
    % Load the data as matrix format for readSHC
    fid = fopen(first_file, 'r');
    if fid == -1
        error('Cannot open file: %s', first_file);
    end
    
    % Skip header lines and read gfc data
    temp_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && startsWith(strtrim(line), 'gfc')
            % Parse gfc line: gfc n m cnm snm sigma_cnm sigma_snm
            parts = str2num(line(4:end)); %#ok<ST2NM>
            if length(parts) >= 4
                temp_data = [temp_data; parts(1:4)]; %#ok<AGROW>
            end
        end
    end
    fclose(fid);
    
    if isempty(temp_data)
        error('No gfc data found in file: %s', first_file);
    end
    
    % Get dimensions
    nmax = max(temp_data(:, 1)); % Maximum degree
    [cnm_temp, snm_temp] = readSHC(temp_data);
    
    fprintf('Maximum degree: %d\n', nmax);
    fprintf('Array dimensions: %d x %d\n', size(cnm_temp, 1), size(cnm_temp, 2));
    
catch ME
    error('Error reading first GRACE file: %s\nError: %s', first_file, ME.message);
end

%% Step 5: Initialize time series arrays
cnm_ts = zeros(nmax + 1, nmax + 1, n_files);
snm_ts = zeros(nmax + 1, nmax + 1, n_files);

%% Step 6: Process each GRACE file
fprintf('Processing %d GRACE files...\n', n_files);

for i = 1:n_files
    if mod(i, 50) == 0 || i == n_files
        fprintf('  Processing file %d/%d (%.1f%%)\n', i, n_files, 100*i/n_files);
    end
    
    current_file = fullfile(grace_dir, gfc_files(i).name);
    current_time = mid_times(i);
    
    try
        % Read GRACE file
        fid = fopen(current_file, 'r');
        temp_data = [];
        
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line) && startsWith(strtrim(line), 'gfc')
                parts = str2num(line(4:end)); %#ok<ST2NM>
                if length(parts) >= 4
                    temp_data = [temp_data; parts(1:4)]; %#ok<AGROW>
                end
            end
        end
        fclose(fid);
        
        if isempty(temp_data)
            warning('No data in file: %s', gfc_files(i).name);
            continue;
        end
        
        % Convert to coefficient matrices
        [cnm, snm] = readSHC(temp_data);
        
        % Ensure proper dimensions
        if size(cnm, 1) ~= nmax + 1 || size(cnm, 2) ~= nmax + 1
            cnm_temp = zeros(nmax + 1, nmax + 1);
            snm_temp = zeros(nmax + 1, nmax + 1);
            
            max_n = min(size(cnm, 1) - 1, nmax);
            cnm_temp(1:max_n+1, 1:max_n+1) = cnm(1:max_n+1, 1:max_n+1);
            snm_temp(1:max_n+1, 1:max_n+1) = snm(1:max_n+1, 1:max_n+1);
            
            cnm = cnm_temp;
            snm = snm_temp;
        end
        
        % Apply C20 replacement
        c20_interp = interp1(c20_time, c20_values, current_time, 'linear', 'extrap');
        cnm(3, 1) = c20_interp; % n=2, m=0 -> index (3,1)
        
        % Add degree-1 coefficients
        if nmax >= 1
            c10_interp = interp1(deg1_time, deg1_c10, current_time, 'linear', 'extrap');
            c11_interp = interp1(deg1_time, deg1_c11, current_time, 'linear', 'extrap');
            s11_interp = interp1(deg1_time, deg1_s11, current_time, 'linear', 'extrap');
            
            % Convert from unnormalized to normalized (multiply by sqrt(3))
            cnm(2, 1) = c10_interp * sqrt(3); % n=1, m=0 -> index (2,1)
            cnm(2, 2) = c11_interp * sqrt(3); % n=1, m=1 -> index (2,2)
            snm(2, 2) = s11_interp * sqrt(3); % n=1, m=1 -> index (2,2)
        end
        
        % Store in time series
        cnm_ts(:, :, i) = cnm;
        snm_ts(:, :, i) = snm;
        
    catch ME
        warning('Error processing file %s: %s', gfc_files(i).name, ME.message);
        % Fill with NaN for this time step
        cnm_ts(:, :, i) = NaN;
        snm_ts(:, :, i) = NaN;
    end
end

%% Step 7: Quality check and summary
n_valid = sum(~any(any(isnan(cnm_ts), 1), 2));
n_invalid = n_files - n_valid;

fprintf('\nProcessing complete!\n');
fprintf('Valid files: %d/%d (%.1f%%)\n', n_valid, n_files, 100*n_valid/n_files);
if n_invalid > 0
    fprintf('Invalid files: %d (filled with NaN)\n', n_invalid);
end

fprintf('Output dimensions:\n');
fprintf('  cnm_ts: %d x %d x %d\n', size(cnm_ts));
fprintf('  snm_ts: %d x %d x %d\n', size(snm_ts));
fprintf('  time_mjd: %d x 1\n', length(time_mjd));

% Display some statistics
fprintf('C20 statistics:\n');
c20_series = squeeze(cnm_ts(3, 1, :));
valid_c20 = c20_series(~isnan(c20_series));
if ~isempty(valid_c20)
    fprintf('  Mean: %.6e\n', mean(valid_c20));
    fprintf('  Std:  %.6e\n', std(valid_c20));
    fprintf('  Range: %.6e to %.6e\n', min(valid_c20), max(valid_c20));
end

fprintf('Degree-1 C10 statistics:\n');
if nmax >= 1
    c10_series = squeeze(cnm_ts(2, 1, :));
    valid_c10 = c10_series(~isnan(c10_series));
    if ~isempty(valid_c10)
        fprintf('  Mean: %.6e\n', mean(valid_c10));
        fprintf('  Std:  %.6e\n', std(valid_c10));
        fprintf('  Range: %.6e to %.6e\n', min(valid_c10), max(valid_c10));
    end
end

end