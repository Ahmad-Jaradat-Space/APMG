function t_quake = select_earthquakes(HistoryMax)
%SELECT_EARTHQUAKES - Selects significant earthquakes by distance and magnitude.
%
% INPUTS:
%   HistoryMean: table with fields 'lat', 'lon', 'mag', 'date'
%   station_lat, station_lon: coordinates of station in degrees
%   max_dist_km: maximum distance in km
%   min_mag: minimum magnitude threshold
%
% OUTPUT:
%   t_quake: column vector of MJD event times (for design matrix)
%



% Extract earthquake info
quake_lat = HistoryMax.latitude;
quake_lon = HistoryMax.longitude;
quake_mag = HistoryMax.mag;



station_lat = 36.1057;
station_lon = 140.0875;
max_dist_km = 500;
min_mag = 7.0;

% % Convert to radians
% lat1 = deg2rad(station_lat);
% lon1 = deg2rad(station_lon);
% lat2 = deg2rad(quake_lat);
% lon2 = deg2rad(quake_lon);
% 
% % Haversine formula
% R = 6371; % Earth radius [km]
% dlat = lat2 - lat1;
% dlon = lon2 - lon1;
% a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
% c = 2*atan2(sqrt(a), sqrt(1-a));
% dist_km = R * c;

dist_deg = distance(station_lat, station_lon, quake_lat, quake_lon);
dist_km = deg2km(dist_deg);


% Selection mask
selected = (dist_km < max_dist_km) & (quake_mag > min_mag);

% Extract dates and convert to MJD
quake_dates = HistoryMax.date(selected);
t_quake = mjuliandate(quake_dates);

end
