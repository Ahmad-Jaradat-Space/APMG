% Debug GPS unit assumptions and conversions
addpath(genpath(pwd));

fprintf('=== GPS UNIT DEBUG ANALYSIS ===\n\n');

% Load GPS coordinates
gps_coords_file = 'data/gps/GPSLatLong.tenv3';
fid = fopen(gps_coords_file, 'r');
if fid == -1
    error('Cannot find GPS coordinates file');
end

station_map = containers.Map();
line_num = 0;
while ~feof(fid)
    line = fgetl(fid);
    line_num = line_num + 1;
    if ischar(line) && ~isempty(line)
        if line_num == 1 && contains(line, 'sta_id')
            continue;
        end
        parts = strsplit(strtrim(line));
        if length(parts) >= 3
            station_id = parts{1};
            lat_val = str2double(parts{2});
            lon_val = str2double(parts{3});
            if ~isnan(lat_val) && ~isnan(lon_val)
                station_map(station_id) = [lat_val, lon_val];
            end
        end
    end
end
fclose(fid);

station_names = keys(station_map);
fprintf('Found %d GPS stations\n\n', length(station_names));

% Test first available station
for i = 1:length(station_names)
    station_file = fullfile('data/gps', sprintf('%s.tenv3', station_names{i}));
    if exist(station_file, 'file')
        fprintf('=== ANALYZING STATION: %s ===\n', station_names{i});
        
        % Read raw file content first
        fprintf('Raw file content (first 10 lines):\n');
        fid = fopen(station_file, 'r');
        for line_count = 1:10
            line = fgetl(fid);
            if ischar(line)
                fprintf('%d: %s\n', line_count, line);
            else
                break;
            end
        end
        fclose(fid);
        
        % Load using existing function
        gps_struct = load_tenv3(station_file);
        
        if ~isempty(gps_struct) && isfield(gps_struct, 'up')
            fprintf('\n=== GPS DATA ANALYSIS ===\n');
            fprintf('Number of observations: %d\n', length(gps_struct.up));
            fprintf('Time range: %.2f to %.2f (MJD)\n', min(gps_struct.t), max(gps_struct.t));
            fprintf('\nVertical (UP) component statistics:\n');
            fprintf('  Raw min: %.6f\n', min(gps_struct.up));
            fprintf('  Raw max: %.6f\n', max(gps_struct.up));
            fprintf('  Raw mean: %.6f\n', mean(gps_struct.up));
            fprintf('  Raw std: %.6f\n', std(gps_struct.up));
            fprintf('  Raw range: %.6f\n', range(gps_struct.up));
            
            % Check if values look like mm, m, or cm
            range_val = range(gps_struct.up);
            mean_abs = mean(abs(gps_struct.up - mean(gps_struct.up)));
            
            fprintf('\n=== UNIT INTERPRETATION ===\n');
            if range_val > 1000
                fprintf('Range > 1000: Likely MILLIMETERS (current assumption)\n');
            elseif range_val > 10 && range_val < 1000
                fprintf('Range 10-1000: Likely MILLIMETERS or large CENTIMETERS\n');
            elseif range_val > 0.1 && range_val < 10
                fprintf('Range 0.1-10: Could be CENTIMETERS or large METERS\n');
            elseif range_val < 0.1
                fprintf('Range < 0.1: Likely METERS (WRONG ASSUMPTION!)\n');
            end
            
            fprintf('\nCurrent processing pipeline:\n');
            fprintf('1. Raw GPS: %.6f to %.6f\n', min(gps_struct.up), max(gps_struct.up));
            
            % Detrend (same as main pipeline)
            poly_result = fitPolynomial(gps_struct.t, gps_struct.up, 1, 0.001);
            detrended = poly_result.v;
            fprintf('2. After detrending: %.6f to %.6f\n', min(detrended), max(detrended));
            
            % Current conversion (assuming mm -> m)
            converted_current = detrended / 1000.0;
            fprintf('3. After /1000 (assume mm->m): %.6f to %.6f\n', min(converted_current), max(converted_current));
            
            % Alternative interpretations
            fprintf('\n=== ALTERNATIVE INTERPRETATIONS ===\n');
            if range_val < 1
                fprintf('If GPS data is already in METERS:\n');
                fprintf('  No conversion needed: %.6f to %.6f m\n', min(detrended), max(detrended));
                fprintf('  In mm this would be: %.2f to %.2f mm\n', min(detrended)*1000, max(detrended)*1000);
            end
            
            if range_val > 1 && range_val < 100
                fprintf('If GPS data is in CENTIMETERS:\n');
                fprintf('  Convert cm->m: %.6f to %.6f m\n', min(detrended)/100, max(detrended)/100);
                fprintf('  In mm this would be: %.2f to %.2f mm\n', min(detrended)*10, max(detrended)*10);
            end
            
        end
        break; % Analyze only first available station
    end
end

fprintf('\n=== RECOMMENDED ACTION ===\n');
fprintf('Check the .tenv3 file format specification!\n');
fprintf('GPS stations typically measure cm-level deformation for drought/seasonal effects.\n');
fprintf('Current assumption may need correction based on actual file format.\n');