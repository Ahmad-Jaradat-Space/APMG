function [lat,lon,time,data,type] = read_oc_netcdf_bsh_u10(filename1)
%READ_OC_NETCDF Summary of this function goes here
%   Detailed explanation goes here
        lat = ncread(filename1,'lat');
        lon = ncread(filename1,'lon');
        time1 = ncread(filename1,'time');
        data = ncread(filename1,'wind_speed');
        type = 1;
       
% convert time in MODEL from days since 1900 to sec1985
    jd1900 = date2jd(1900);
    jd1985  = date2jd(1985);
    time = (double(time1) - jd1985 + jd1900)*86400;
end

