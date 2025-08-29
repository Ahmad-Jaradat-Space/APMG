function [cnm_static, snm_static] = computeStaticReferenceField(grace_dir, c20_file, deg1_file, reference_years, nmax)
% Use BA01 files with specified nmax to match main processing
if nargin < 5
    nmax = 60;  % Default to degree 60 for BA01 files
end
gfc_files = dir(fullfile(grace_dir, '*BA01*.gfc'));
reference_files = {};
reference_dates = [];
for i = 1:length(gfc_files)
    filename = gfc_files(i).name;
    pattern = 'GSM-2_(\d{7})-(\d{7})_';
    tokens = regexp(filename, pattern, 'tokens');
    if ~isempty(tokens)
        start_date = str2double(tokens{1}{1});
        end_date = str2double(tokens{1}{2});
        start_year = floor(start_date / 1000);
        end_year = floor(end_date / 1000);
        if any(start_year == reference_years) || any(end_year == reference_years)
            reference_files{end+1} = filename;
            start_doy = mod(start_date, 1000);
            end_doy = mod(end_date, 1000);
            mid_year = (start_year + end_year) / 2;
            mid_decimal_year = mid_year + (start_doy + end_doy) / (2 * 365.25);
            reference_dates(end+1) = mid_decimal_year;
        end
    end
end
n_files = length(reference_files);
[reference_dates, sort_idx] = sort(reference_dates);
reference_files = reference_files(sort_idx);
fid = fopen(c20_file, 'r');
c20_data = [];
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~startsWith(strtrim(line), '#') && ~isempty(strtrim(line))
        nums = str2num(line);
        if ~isempty(nums) && length(nums) >= 2
            c20_data = [c20_data; nums(1:2)];
        end
    end
end
fclose(fid);
c20_time = c20_data(:, 1);
c20_values = c20_data(:, 2);
first_file = fullfile(grace_dir, reference_files{1});
fid = fopen(first_file, 'r');
temp_data = [];
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && startsWith(strtrim(line), 'gfc')
        parts = str2num(line(4:end));
        if length(parts) >= 4
            temp_data = [temp_data; parts(1:4)];
        end
    end
end
fclose(fid);
cnm_sum = zeros(nmax + 1, nmax + 1);
snm_sum = zeros(nmax + 1, nmax + 1);
valid_count = 0;
for i = 1:n_files
    current_file = fullfile(grace_dir, reference_files{i});
    current_time = reference_dates(i);
    fid = fopen(current_file, 'r');
    temp_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && startsWith(strtrim(line), 'gfc')
            parts = str2num(line(4:end));
            if length(parts) >= 4
                temp_data = [temp_data; parts(1:4)];
            end
        end
    end
    fclose(fid);
    [cnm, snm] = readSHC(temp_data, nmax);
    c20_interp = interp1(c20_time, c20_values, current_time, 'linear', 'extrap');
    cnm(3, 1) = c20_interp;
    cnm_sum = cnm_sum + cnm;
    snm_sum = snm_sum + snm;
    valid_count = valid_count + 1;
end
cnm_static = cnm_sum / valid_count;
snm_static = snm_sum / valid_count;
end