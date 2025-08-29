function [cnm_ts, snm_ts, time_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file)
% Use only BB01 files (degree 96, higher resolution)
gfc_files = dir(fullfile(grace_dir, '*BB01*.gfc'));

n_files = length(gfc_files);
file_dates = zeros(n_files, 2);
for i = 1:n_files
    filename = gfc_files(i).name;
    pattern = 'GSM-2_(\d{7})-(\d{7})_';
    tokens = regexp(filename, pattern, 'tokens');
    if ~isempty(tokens)
        start_date = str2double(tokens{1}{1});
        end_date = str2double(tokens{1}{2});
        start_year = floor(start_date / 1000);
        start_doy = mod(start_date, 1000);
        end_year = floor(end_date / 1000);
        end_doy = mod(end_date, 1000);
        file_dates(i, 1) = start_year + (start_doy - 1) / 365.25;
        file_dates(i, 2) = end_year + (end_doy - 1) / 365.25;
    end
end
mid_times = mean(file_dates, 2);
time_mjd = decyear2mjd(mid_times);
[time_mjd, sort_idx] = sort(time_mjd);
gfc_files = gfc_files(sort_idx);

% Load C20 replacement (normalized)
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

% Load degree-1 normalized coefficients from file with explicit n,m tags
fid = fopen(deg1_file, 'r');
raw_rows = [];
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    if isempty(strtrim(line)) || startsWith(strtrim(line), 'TITLE') || startsWith(strtrim(line), 'UPDATE') || startsWith(strtrim(line), 'AUTHOR') || startsWith(strtrim(line), 'DESCRIPTION') || startsWith(strtrim(line), 'We have') || startsWith(strtrim(line), 'For ocean') || startsWith(strtrim(line), 'Data format') || startsWith(strtrim(line), 'REFERENCE') || startsWith(strtrim(line), 'SPECIAL') || startsWith(strtrim(line), '''')
        continue;
    end
    nums = str2num(line); %#ok<ST2NM>
    % Expect format: YYYYMM  n  m   C   S  ... start end
    if ~isempty(nums) && length(nums) >= 5
        raw_rows = [raw_rows; nums(1:5)]; %#ok<AGROW>
    end
end
fclose(fid);

% Split by (n,m)
mask10 = raw_rows(:,2)==1 & raw_rows(:,3)==0;
mask11 = raw_rows(:,2)==1 & raw_rows(:,3)==1;
Y10 = raw_rows(mask10,1); C10 = raw_rows(mask10,4); S10 = raw_rows(mask10,5);
Y11 = raw_rows(mask11,1); C11 = raw_rows(mask11,4); S11 = raw_rows(mask11,5);

% Convert YYYYMM to decimal year mid-month
convYM = @(ym) (floor(ym/100) + (mod(ym,100)-0.5)/12);
t10 = arrayfun(convYM, Y10);
t11 = arrayfun(convYM, Y11);

% Read first file to get nmax
first_file = fullfile(grace_dir, gfc_files(1).name);
fid = fopen(first_file, 'r');
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
cnm_ts = zeros(nmax + 1, nmax + 1, n_files);
snm_ts = zeros(nmax + 1, nmax + 1, n_files);

for i = 1:n_files
    current_file = fullfile(grace_dir, gfc_files(i).name);
    current_time = mjd2decyear(time_mjd(i));
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
    [cnm, snm] = readSHC(temp_data);

    % C20 replacement (normalized)
    c20_interp = interp1(c20_time, c20_values, current_time, 'linear', 'extrap');
    cnm(3, 1) = c20_interp;

    % Degree-1 (normalized, no sqrt(3))
    if nmax >= 1
        if ~isempty(t10)
            c10_interp = interp1(t10, C10, current_time, 'linear', 'extrap');
            % S10 often ~0 in product; we don't have a slot for S10 (m=0 sine term is zero by definition)
            cnm(2, 1) = c10_interp;
        end
        if ~isempty(t11)
            c11_interp = interp1(t11, C11, current_time, 'linear', 'extrap');
            s11_interp = interp1(t11, S11, current_time, 'linear', 'extrap');
            cnm(2, 2) = c11_interp;
            snm(2, 2) = s11_interp;
        end
    end

    cnm_ts(:, :, i) = cnm;
    snm_ts(:, :, i) = snm;
end
end