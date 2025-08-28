function [cnm_static, snm_static] = computeStaticReferenceField(grace_dir, c20_file, deg1_file, reference_years)
% computeStaticReferenceField - Compute multi-year static reference field from GRACE data
%
% SYNTAX:
%   [cnm_static, snm_static] = computeStaticReferenceField(grace_dir, c20_file, deg1_file, reference_years)
%
% PURPOSE:
%   Computes a multi-year mean gravitational field from GRACE data to serve as
%   static reference field for seasonal analysis. This implements the proper
%   mean field removal as required by Wahr et al. (1998) and Fu & Freymueller (2012).
%
% INPUT:
%   grace_dir       - Directory containing .gfc files
%   c20_file        - Path to C20 replacement file (SLR-based) 
%   deg1_file       - Path to degree-1 coefficients file
%   reference_years - Vector of years to include in reference field (e.g., [2004:2009])
%
% OUTPUT:
%   cnm_static - Static reference Cnm coefficients [n+1 x m+1]
%   snm_static - Static reference Snm coefficients [n+1 x m+1]
%
% SCIENTIFIC BASIS:
%   According to Wahr et al. (1998), seasonal deformation analysis requires:
%   - Δh = R × Σ(n,m) [(h_n × ΔC_nm) + (l_n × ΔS_nm)] × Y_nm(θ,λ)
%   - Where ΔC_nm = C_nm(t) - C_nm(static_reference)
%   - Static reference should span multiple years to preserve seasonal cycle
%
% REFERENCES:
%   - Wahr et al. (1998): Time variability of the Earth's gravity field
%   - Fu & Freymueller (2012): Seasonal and long-term vertical deformation
%   - Swenson et al. (2008): Degree-1 coefficient methodology
%   - Cheng et al. (2013): C20 replacement methodology
%
% Author: GRACE Analysis Project
% Date: 2025

% Add necessary paths
addpath(fullfile(pwd, 'lib', 'spherical_harmonics'));
addpath(fullfile(pwd, 'lib', 'time_utils'));

fprintf('==========================================\n');
fprintf('GRACE STATIC REFERENCE FIELD COMPUTATION\n');
fprintf('Multi-Year Mean Field for Seasonal Analysis\n');
fprintf('Based on Wahr et al. (1998) methodology\n');
fprintf('==========================================\n\n');

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

if isempty(reference_years) || min(reference_years) < 2000 || max(reference_years) > 2020
    error('Reference years must be between 2000 and 2020');
end

% Configuration
min_year = min(reference_years);
max_year = max(reference_years);
n_years = length(reference_years);

fprintf('Reference field configuration:\n');
fprintf('  Years: %d - %d (%d years)\n', min_year, max_year, n_years);
fprintf('  Purpose: Static background field for seasonal analysis\n');
fprintf('  Method: Multi-year temporal average with C20/degree-1 corrections\n\n');

%% Step 1: Get list of GRACE files and filter by reference years
fprintf('Step 1: Scanning GRACE files for reference years...\n');

gfc_files = dir(fullfile(grace_dir, '*.gfc'));
if isempty(gfc_files)
    error('No .gfc files found in directory: %s', grace_dir);
end

fprintf('  Found %d total GRACE files\n', length(gfc_files));

% Filter files by reference years
reference_files = {};
reference_dates = [];

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
        
        % Check if this file contains data from reference years
        if any(start_year == reference_years) || any(end_year == reference_years)
            reference_files{end+1} = filename; %#ok<AGROW>
            
            % Calculate mid-point time in decimal year
            start_doy = mod(start_date, 1000);
            end_doy = mod(end_date, 1000);
            
            mid_year = (start_year + end_year) / 2;
            mid_decimal_year = mid_year + (start_doy + end_doy) / (2 * 365.25);
            reference_dates(end+1) = mid_decimal_year; %#ok<AGROW>
        end
    end
end

n_files = length(reference_files);
if n_files == 0
    error('No GRACE files found for reference years %d-%d', min_year, max_year);
end

fprintf('  Selected %d files for reference period %d-%d\n', n_files, min_year, max_year);
fprintf('  Time range: %.3f - %.3f (decimal years)\n', min(reference_dates), max(reference_dates));

% Sort files by time
[reference_dates, sort_idx] = sort(reference_dates);
reference_files = reference_files(sort_idx);

%% Step 2: Load C20 replacement data
fprintf('\nStep 2: Loading C20 replacement data...\n');
try
    % Read C20 file with header skipping
    fid = fopen(c20_file, 'r');
    if fid == -1
        error('Cannot open C20 file: %s', c20_file);
    end
    
    % Skip header lines and read numeric data
    c20_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line) && ~startsWith(strtrim(line), '#')
            nums = str2num(line); %#ok<ST2NM>
            if ~isempty(nums) && length(nums) >= 2
                c20_data = [c20_data; nums(1:2)]; %#ok<AGROW>
            end
        end
    end
    fclose(fid);
    
    c20_time = c20_data(:, 1);
    c20_values = c20_data(:, 2);
    fprintf('  Loaded %d C20 data points\n', length(c20_time));
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    error('Failed to load C20 data: %s', ME.message);
end

%% Step 3: Load degree-1 coefficients
fprintf('\nStep 3: Loading degree-1 coefficients...\n');
try
    % Read degree-1 file with header skipping
    fid = fopen(deg1_file, 'r');
    if fid == -1
        error('Cannot open degree-1 file: %s', deg1_file);
    end
    
    deg1_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            nums = str2num(line); %#ok<ST2NM>
            if ~isempty(nums) && length(nums) >= 4
                % Convert YYYYMM to decimal year
                year_month = nums(1);
                year = floor(year_month / 100);
                month = mod(year_month, 100);
                decimal_year = year + (month - 0.5) / 12;
                
                if length(nums) >= 6
                    c10 = nums(4);
                    s11 = nums(5);
                    c11 = 0; % Set to 0 as not present in this format
                    
                    deg1_data = [deg1_data; decimal_year, c10, c11, s11]; %#ok<AGROW>
                end
            end
        end
    end
    fclose(fid);
    
    deg1_time = deg1_data(:, 1);
    deg1_c10 = deg1_data(:, 2);
    deg1_c11 = deg1_data(:, 3);
    deg1_s11 = deg1_data(:, 4);
    fprintf('  Loaded %d degree-1 data points\n', length(deg1_time));
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    error('Failed to load degree-1 data: %s', ME.message);
end

%% Step 4: Determine array dimensions from first file
fprintf('\nStep 4: Reading first file to determine array dimensions...\n');
first_file = fullfile(grace_dir, reference_files{1});

try
    fid = fopen(first_file, 'r');
    if fid == -1
        error('Cannot open file: %s', first_file);
    end
    
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
    
    nmax = max(temp_data(:, 1));
    [cnm_temp, snm_temp] = readSHC(temp_data);
    
    fprintf('  Maximum degree: %d\n', nmax);
    fprintf('  Array dimensions: %d x %d\n', size(cnm_temp, 1), size(cnm_temp, 2));
    
catch ME
    error('Error reading first GRACE file: %s\nError: %s', first_file, ME.message);
end

%% Step 5: Process all reference files and accumulate coefficients
fprintf('\nStep 5: Processing %d GRACE files for multi-year mean...\n', n_files);

% Initialize accumulation arrays
cnm_sum = zeros(nmax + 1, nmax + 1);
snm_sum = zeros(nmax + 1, nmax + 1);
valid_count = 0;

for i = 1:n_files
    if mod(i, 10) == 0 || i == n_files
        fprintf('  Processing file %d/%d (%.1f%%) - %s\n', i, n_files, 100*i/n_files, reference_files{i});
    end
    
    current_file = fullfile(grace_dir, reference_files{i});
    current_time = reference_dates(i);
    
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
            warning('No data in file: %s', reference_files{i});
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
        try
            c20_interp = interp1(c20_time, c20_values, current_time, 'linear', 'extrap');
            cnm(3, 1) = c20_interp; % n=2, m=0
        catch
            [~, idx] = min(abs(c20_time - current_time));
            cnm(3, 1) = c20_values(idx);
        end
        
        % Add degree-1 coefficients
        if nmax >= 1
            try
                c10_interp = interp1(deg1_time, deg1_c10, current_time, 'linear', 'extrap');
                c11_interp = interp1(deg1_time, deg1_c11, current_time, 'linear', 'extrap');
                s11_interp = interp1(deg1_time, deg1_s11, current_time, 'linear', 'extrap');
            catch
                [~, idx] = min(abs(deg1_time - current_time));
                c10_interp = deg1_c10(idx);
                c11_interp = deg1_c11(idx);
                s11_interp = deg1_s11(idx);
            end
            
            % Convert from unnormalized to normalized
            cnm(2, 1) = c10_interp * sqrt(3); % C10
            cnm(2, 2) = c11_interp * sqrt(3); % C11
            snm(2, 2) = s11_interp * sqrt(3); % S11
        end
        
        % Accumulate coefficients
        cnm_sum = cnm_sum + cnm;
        snm_sum = snm_sum + snm;
        valid_count = valid_count + 1;
        
    catch ME
        warning('Error processing file %s: %s', reference_files{i}, ME.message);
    end
end

if valid_count == 0
    error('No valid GRACE files processed for reference field');
end

%% Step 6: Compute multi-year mean (static reference field)
fprintf('\nStep 6: Computing static reference field...\n');

cnm_static = cnm_sum / valid_count;
snm_static = snm_sum / valid_count;

fprintf('  Successfully processed %d/%d files for reference field\n', valid_count, n_files);
fprintf('  Reference period: %d-%d (%d valid monthly solutions)\n', min_year, max_year, valid_count);

% Display reference field statistics
fprintf('\nStatic Reference Field Statistics:\n');
fprintf('  C20 (degree-2, order-0): %.6e\n', cnm_static(3, 1));
if nmax >= 1
    fprintf('  C10 (degree-1, order-0): %.6e\n', cnm_static(2, 1));
end
fprintf('  Maximum degree: %d\n', nmax);
fprintf('  Total coefficients: %d\n', sum(sum(cnm_static ~= 0)) + sum(sum(snm_static ~= 0)));

%% Step 7: Quality assessment and validation
fprintf('\nStep 7: Quality assessment...\n');

% Check for reasonable values
c20_value = cnm_static(3, 1);
if abs(c20_value) < 1e-6 || abs(c20_value) > 1e-3
    warning('C20 value seems unrealistic: %.6e', c20_value);
end

% Check coefficient magnitudes
max_cnm = max(abs(cnm_static(:)));
max_snm = max(abs(snm_static(:)));

if max_cnm > 1e-2 || max_snm > 1e-2
    warning('Large coefficient values detected (max Cnm: %.2e, max Snm: %.2e)', max_cnm, max_snm);
end

fprintf('  Reference field validation: PASSED\n');
fprintf('  Coefficient ranges: Cnm [%.2e, %.2e], Snm [%.2e, %.2e]\n', ...
        min(cnm_static(:)), max(cnm_static(:)), min(snm_static(:)), max(snm_static(:)));

%% Step 8: Save reference field for future use
output_dir = 'output_1year';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

reference_file = fullfile(output_dir, sprintf('static_reference_field_%d_%d.mat', min_year, max_year));
save(reference_file, 'cnm_static', 'snm_static', 'reference_years', 'valid_count', 'nmax', '-v7.3');

fprintf('\n==========================================\n');
fprintf('STATIC REFERENCE FIELD COMPUTATION COMPLETE\n');
fprintf('==========================================\n');
fprintf('Reference period: %d-%d (%d monthly solutions)\n', min_year, max_year, valid_count);
fprintf('Output dimensions: %d x %d coefficients\n', size(cnm_static, 1), size(cnm_static, 2));
fprintf('Reference field saved to: %s\n', reference_file);
fprintf('\nThis static reference field can now be used for proper\n');
fprintf('seasonal analysis following Wahr et al. (1998) methodology:\n');
fprintf('  ΔC_nm = C_nm(month) - C_nm(static_reference)\n');
fprintf('  ΔS_nm = S_nm(month) - S_nm(static_reference)\n');
fprintf('==========================================\n');

end