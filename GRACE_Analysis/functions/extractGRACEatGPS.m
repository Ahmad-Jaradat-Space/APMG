function grace_at_gps = extractGRACEatGPS(u_vertical_grid, lat_grid, lon_grid, lat_gps, lon_gps)
if ndims(u_vertical_grid) == 3
    [nlat, nlon, ntime] = size(u_vertical_grid);
else
    [nlat, nlon] = size(u_vertical_grid);
    ntime = 1;
end
nstations = length(lat_gps);
grace_at_gps = zeros(nstations, ntime);
for t = 1:ntime
    if ndims(u_vertical_grid) == 3
        current_deformation = u_vertical_grid(:, :, t);
    else
        current_deformation = u_vertical_grid;
    end
    for s = 1:nstations
        interpolated_value = interp2(lon_grid, lat_grid, current_deformation, lon_gps(s), lat_gps(s), 'linear');
        grace_at_gps(s, t) = interpolated_value;
    end
end
end