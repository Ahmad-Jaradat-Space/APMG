function dist = get_dist2station(lat_data,lon_data,lat_station,lon_station)
%GET_DIST2STATION Summary of this function goes here
%   Detailed explanation goes here
R = 6371.000785; % mean Earthradius [km}
dist_north = (lat_data - lat_station)*pi/180*R;

R_lat = cosd(lat_data)*R;
dist_west = (lon_data - lon_station)*pi/180.*R_lat;

dist = sqrt(dist_north.^2 + dist_west.^2);

end

