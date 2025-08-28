function [dist2coast] = betterDist2Coast_lola(lon2,lat2)

%% Autor: Shaun Richardt
% Input: Matrix  where Latitude is in (:,2), Longitude (:,3) and dist2coast (:,24)
% Output: same Matrix but with updated dist2coast (:,24) (metres)
% %Filepath from the Grid @line 11
% data = readNetCDF2('path to the grid.nc');

%%
%Filepath from the Grid
%data = readNetCDF2([pwd '/../../../../data/coastline/dist_to_GSHHG_v2.3.4_1m.nc']);
data = readNetCDF2('/Users/fenoglio/sciebo/shaun/data/coastline/dist_to_GSHHG_v2.3.4_1m.nc');

% lat2 = altimeter(:,2);
% lon2 = altimeter(:,3);

%interp2(X,Y,Z,XI,YI);
distcoast = interp2(data.lat, data.lon, data.z, lat2, lon2);

dist2coast = distcoast * -1000; % in m

%altimeter(:,24) = distcoast;

end


