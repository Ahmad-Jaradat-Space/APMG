function [lat_gps, lon_gps, gps_data, time_mjd_gps, station_names] = readGPS(gps_coords_file, gps_data_dir)

% Read coordinates file
fid = fopen(gps_coords_file, 'r');
if fid == -1
    return;
end
fgetl(fid); % skip header
data = textscan(fid, '%s %f %f', 'CollectOutput', false);
fclose(fid);

names = data{1};
lats = data{2};  
lons = data{3};

% Load GPS files
lat_gps = [];
lon_gps = [];
gps_data = {};
time_mjd_gps = {};
station_names = {};

for i = 1:length(names)
    station_file = fullfile(gps_data_dir, sprintf('%s.tenv3', names{i}));
    if exist(station_file, 'file')
        gps_struct = load_tenv3(station_file);
        elevation = gps_struct.up;
        
        if strcmp(names{i}, 'P056')
            time_decyear = mjd2decyear(gps_struct.t);
            poly_result = fitPiecewiseLinear(time_decyear, elevation, 2014.073);
            elevation_detrended = poly_result.v;
        else
            poly_result = fitPolynomial(gps_struct.t, elevation, 1, 0.001);
            elevation_detrended = poly_result.v;
        end
        
        lat_gps(end+1) = lats(i);
        lon_gps(end+1) = lons(i);
        gps_data{end+1} = elevation_detrended;
        time_mjd_gps{end+1} = gps_struct.t;
        station_names{end+1} = names{i};
    end
end

lat_gps = lat_gps(:);
lon_gps = lon_gps(:);

end