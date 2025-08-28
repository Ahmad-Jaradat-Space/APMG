function [lat,lon,time,data,type] = read_oc_netcdf_bsh_tid(filename)
%READ_OC_NETCDF Summary of this function goes here
%   Detailed explanation goes here
        lat = ncread(filename,'ytno');
        lon = ncread(filename,'xtno');
        time1 = ncread(filename,'time');
        data = ncread(filename,'ztideno');
        type = 4;
    
      
% convert time in MODEL from days since 1900 to sec1985
    jd1900 = date2jd(1900);
    jd1985  = date2jd(1985);
    daydiff = jd1900 - jd1985;
    time = (double(time1) - jd1985 + jd1900)*86400;
    
end

