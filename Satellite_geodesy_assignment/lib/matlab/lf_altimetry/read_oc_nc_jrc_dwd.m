function [lat,lon,time,data,type] = read_oc_nc_jrc_dwd(filename)
%READ_OC_NETCDF Summary of this function goes here
%   Detailed explanation goes here
        lat = ncread(filename,'LAT');
        lon = ncread(filename,'LON');
        time1 = ncread(filename,'TIME');
        data = ncread(filename,'HA');
        type = 3;
    
      
% convert time in MODEL from sec since start to sec1985
% start at 28 Nov 2013 00:00 Alessandro 13.11.2014
    von=date2mjd(2013,11,28,00,00,00);
    secstart = mjd_to_sec1985(von);
    time = (double(time1) + secstart);
    
end

