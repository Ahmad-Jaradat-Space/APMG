function [cnm_ts, snm_ts, time_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file)
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
fid = fopen(deg1_file, 'r');
deg1_data = [];
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(strtrim(line))
        nums = str2num(line);
        if ~isempty(nums) && length(nums) >= 5
            year_month = nums(1);
            year = floor(year_month / 100);
            month = mod(year_month, 100);
            decimal_year = year + (month - 0.5) / 12;
            c10 = nums(4);
            s11 = nums(5);
            deg1_data = [deg1_data; decimal_year, c10, 0, s11];
        end
    end
end
fclose(fid);

% Extract degree-1 data arrays for interpolation (Swenson et al. 2008)
deg1_time = deg1_data(:, 1); % Decimal year
deg1_c10 = deg1_data(:, 2);  % C10 (unnormalized)
deg1_c11 = deg1_data(:, 3);  % C11 (unnormalized) - set to 0
deg1_s11 = deg1_data(:, 4);  % S11 (unnormalized)

first_file = fullfile(grace_dir, gfc_files(1).name);
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
            parts = str2num(line(4:end));
            if length(parts) >= 4
                temp_data = [temp_data; parts(1:4)];
            end
        end
    end
    fclose(fid);
    [cnm, snm] = readSHC(temp_data);
    
    % Apply C20 replacement (Cheng et al. 2013)
    c20_interp = interp1(c20_time, c20_values, current_time, 'linear', 'extrap');
    cnm(3, 1) = c20_interp;
    
    % Add degree-1 coefficients (Swenson et al. 2008) - CRITICAL for seasonal analysis
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
    
    cnm_ts(:, :, i) = cnm;
    snm_ts(:, :, i) = snm;
end
end