function [cnm_ts, snm_ts, time_mjd] = processGRACEfiles_1year(grace_dir, c20_file, deg1_file, target_year)
% processGRACEfiles_1year - Process GRACE files for a specific year with corrections
%
% SYNTAX:
%   [cnm_ts, snm_ts, time_mjd] = processGRACEfiles_1year(grace_dir, c20_file, deg1_file, target_year)
%
% INPUT:
%   grace_dir   - Directory containing .gfc files
%   c20_file    - Path to C20 replacement file (SLR-based)
%   deg1_file   - Path to degree-1 coefficients file
%   target_year - Target year for filtering (e.g., 2005)
%
% OUTPUT:
%   cnm_ts    - Time series of Cnm coefficients [n+1 x m+1 x n_months]
%   snm_ts    - Time series of Snm coefficients [n+1 x m+1 x n_months]
%   time_mjd  - Time stamps in Modified Julian Date [n_months x 1]
%
% ENHANCED FEATURES:
%   - Year-based filtering for focused analysis
%   - Improved file name parsing
%   - Enhanced error handling and validation
%
% REFERENCES:
%   - Swenson et al. (2008) for degree-1 coefficients
%   - Cheng et al. (2013) for C20 replacement
%   - Wahr et al. (1998) for GRACE processing methodology
%
% Author: GRACE Analysis Project - Enhanced for 1-year analysis
% Date: 2025

% Add path to spherical harmonics functions
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));

% Input validation
if nargin < 4
    error('All 4 input arguments are required');
end

if ~exist(grace_dir, 'dir')
    error('GRACE directory does not exist: %s', grace_dir);
end

if ~exist(c20_file, 'file')
    error('C20 file does not exist: %s', c20_file);
end

if ~exist(deg1_file, 'file')
    error('Degree-1 file does not exist: %s', deg1_file);
end

if target_year < 2000 || target_year > 2020
    error('Target year must be between 2000 and 2020');
end

fprintf('Processing GRACE files for year %d\n', target_year);
fprintf('GRACE directory: %s\n', grace_dir);
fprintf('C20 corrections: %s\n', c20_file);
fprintf('Degree-1 coefficients: %s\n', deg1_file);

%% Step 1: Get list of GRACE files and filter by year
gfc_files = dir(fullfile(grace_dir, '*.gfc'));
if isempty(gfc_files)
    error('No .gfc files found in directory: %s', grace_dir);
end

fprintf('Found %d GRACE files total\n', length(gfc_files));

% Filter files by target year
target_files = {};
target_dates = [];

for i = 1:length(gfc_files)
    filename = gfc_files(i).name;
    
    % Extract date range from filename using regex
    % Format: kfilter_DDK3_GSM-2_YYYYDDD-YYYYDDD_GRAC_UTCSR_BA01_0600.gfc
    pattern = 'GSM-2_(\d{7})-(\d{7})_';
    tokens = regexp(filename, pattern, 'tokens');
    
    if ~isempty(tokens)
        start_date = str2double(tokens{1}{1});
        end_date = str2double(tokens{1}{2});
        
        % Extract years
        start_year = floor(start_date / 1000);
        end_year = floor(end_date / 1000);
        
        % Check if this file contains data from target year
        if start_year == target_year || end_year == target_year
            target_files{end+1} = filename; %#ok<AGROW>
            
            % Calculate mid-point time in decimal year
            start_doy = mod(start_date, 1000);
            end_doy = mod(end_date, 1000);
            
            mid_decimal_year = target_year + (start_doy + end_doy) / (2 * 365.25);
            target_dates(end+1) = mid_decimal_year; %#ok<AGROW>
        end
    end
end

n_files = length(target_files);
if n_files == 0
    error('No GRACE files found for year %d', target_year);
end

fprintf('Found %d GRACE files for year %d\n', n_files, target_year);

% Convert decimal years to MJD
time_mjd = zeros(n_files, 1);
for i = 1:n_files
    try
        time_mjd(i) = decyear2mjd(target_dates(i));
    catch
        % Fallback calculation if function not available
        time_mjd(i) = (target_dates(i) - 1858.0) * 365.25; % Approximate MJD
    end
end

% Sort files by time
[time_mjd, sort_idx] = sort(time_mjd);
target_files = target_files(sort_idx);
target_dates = target_dates(sort_idx);

fprintf('Time range: %.3f - %.3f (decimal years)\n', min(target_dates), max(target_dates));

%% Step 2: Load C20 replacement data
fprintf('Loading C20 replacement data...\n');
try
    % Read C20 file with header skipping
    fid = fopen(c20_file, 'r');
    if fid == -1
        error('Cannot open C20 file: %s', c20_file);
    end
    
    % Skip header lines (lines starting with # or containing text)
    c20_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line) && ~startsWith(strtrim(line), '#')
            % Try to parse as numeric data
            nums = str2num(line); %#ok<ST2NM>
            if ~isempty(nums) && length(nums) >= 2
                c20_data = [c20_data; nums(1:2)]; %#ok<AGROW>
            end
        end
    end
    fclose(fid);
    
    if isempty(c20_data)
        error('No numeric C20 data found in file');
    end
    
    c20_time = c20_data(:, 1); % Decimal year
    c20_values = c20_data(:, 2); % Normalized C20 values
    fprintf('  Loaded %d C20 data points\n', length(c20_time));
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    error('Failed to load C20 data: %s', ME.message);
end

%% Step 3: Load degree-1 coefficients
fprintf('Loading degree-1 coefficients...\n');
try
    % Read degree-1 file with header skipping
    fid = fopen(deg1_file, 'r');
    if fid == -1
        error('Cannot open degree-1 file: %s', deg1_file);
    end
    
    % Skip header lines and read numeric data
    deg1_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            % Try to parse as numeric data (should have YYYYMM format)
            nums = str2num(line); %#ok<ST2NM>
            if ~isempty(nums) && length(nums) >= 4
                % Convert YYYYMM to decimal year
                year_month = nums(1);
                year = floor(year_month / 100);
                month = mod(year_month, 100);
                decimal_year = year + (month - 0.5) / 12; % Mid-month
                
                % Extract coefficients (columns 4, 5 seem to be C10, S11 based on format)
                if length(nums) >= 6
                    c10 = nums(4); % C10 coefficient
                    s11 = nums(5); % S11 coefficient  
                    c11 = 0; % Not present in this format, set to 0
                    
                    deg1_data = [deg1_data; decimal_year, c10, c11, s11]; %#ok<AGROW>
                end
            end
        end
    end
    fclose(fid);
    
    if isempty(deg1_data)
        error('No numeric degree-1 data found in file');
    end
    
    deg1_time = deg1_data(:, 1); % Decimal year
    deg1_c10 = deg1_data(:, 2); % C10 (unnormalized)
    deg1_c11 = deg1_data(:, 3); % C11 (unnormalized) - set to 0
    deg1_s11 = deg1_data(:, 4); % S11 (unnormalized)
    fprintf('  Loaded %d degree-1 data points\n', length(deg1_time));
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    error('Failed to load degree-1 data: %s', ME.message);
end

%% Step 4: Read first file to get dimensions
fprintf('Reading first file to determine array dimensions...\n');
first_file = fullfile(grace_dir, target_files{1});

try
    % Read using existing readSHC function
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
fprintf('Processing %d GRACE files for %d...\n', n_files, target_year);

for i = 1:n_files
    if mod(i, 5) == 0 || i == n_files
        fprintf('  Processing file %d/%d (%.1f%%) - %s\n', i, n_files, 100*i/n_files, target_files{i});
    end
    
    current_file = fullfile(grace_dir, target_files{i});
    current_time = target_dates(i);
    
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
            warning('No data in file: %s', target_files{i});
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
        
        % Apply C20 replacement (with error handling)
        try
            c20_interp = interp1(c20_time, c20_values, current_time, 'linear', 'extrap');
            cnm(3, 1) = c20_interp; % n=2, m=0 -> index (3,1)
        catch
            % Use nearest neighbor if interpolation fails
            [~, idx] = min(abs(c20_time - current_time));
            cnm(3, 1) = c20_values(idx);
        end
        
        % Add degree-1 coefficients (with error handling)
        if nmax >= 1
            try
                c10_interp = interp1(deg1_time, deg1_c10, current_time, 'linear', 'extrap');
                c11_interp = interp1(deg1_time, deg1_c11, current_time, 'linear', 'extrap');
                s11_interp = interp1(deg1_time, deg1_s11, current_time, 'linear', 'extrap');
            catch
                % Use nearest neighbor if interpolation fails
                [~, idx] = min(abs(deg1_time - current_time));
                c10_interp = deg1_c10(idx);
                c11_interp = deg1_c11(idx);
                s11_interp = deg1_s11(idx);
            end
            
            % Convert from unnormalized to normalized (multiply by sqrt(3))
            cnm(2, 1) = c10_interp * sqrt(3); % n=1, m=0 -> index (2,1)
            cnm(2, 2) = c11_interp * sqrt(3); % n=1, m=1 -> index (2,2)
            snm(2, 2) = s11_interp * sqrt(3); % n=1, m=1 -> index (2,2)
        end
        
        % Store in time series
        cnm_ts(:, :, i) = cnm;
        snm_ts(:, :, i) = snm;
        
    catch ME
        warning('Error processing file %s: %s', target_files{i}, ME.message);
        % Fill with NaN for this time step
        cnm_ts(:, :, i) = NaN;
        snm_ts(:, :, i) = NaN;
    end
end

%% Step 7: Quality check and summary
n_valid = sum(~any(any(isnan(cnm_ts), 1), 2));
n_invalid = n_files - n_valid;

fprintf('\nProcessing complete for year %d!\n', target_year);
fprintf('Valid files: %d/%d (%.1f%%)\n', n_valid, n_files, 100*n_valid/n_files);
if n_invalid > 0
    fprintf('Invalid files: %d (filled with NaN)\n', n_invalid);
end

fprintf('Output dimensions:\n');
fprintf('  cnm_ts: %d x %d x %d\n', size(cnm_ts));
fprintf('  snm_ts: %d x %d x %d\n', size(snm_ts));
fprintf('  time_mjd: %d x 1\n', length(time_mjd));

% Display some statistics
fprintf('C20 statistics for %d:\n', target_year);
c20_series = squeeze(cnm_ts(3, 1, :));
valid_c20 = c20_series(~isnan(c20_series));
if ~isempty(valid_c20)
    fprintf('  Mean: %.6e\n', mean(valid_c20));
    fprintf('  Std:  %.6e\n', std(valid_c20));
    fprintf('  Range: %.6e to %.6e\n', min(valid_c20), max(valid_c20));
end

fprintf('Degree-1 C10 statistics for %d:\n', target_year);
if nmax >= 1
    c10_series = squeeze(cnm_ts(2, 1, :));
    valid_c10 = c10_series(~isnan(c10_series));
    if ~isempty(valid_c10)
        fprintf('  Mean: %.6e\n', mean(valid_c10));
        fprintf('  Std:  %.6e\n', std(valid_c10));
        fprintf('  Range: %.6e to %.6e\n', min(valid_c10), max(valid_c10));
    end
end

fprintf('Time coverage for %d: MJD %.1f to %.1f\n', target_year, min(time_mjd), max(time_mjd));

end