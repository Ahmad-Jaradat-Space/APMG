function [tsec85] = day50_to_sec85(tday50)
%READ_OC_NETCDF Summary of this function goes here
% convert time from days since 1950 to sec1985
    jd1950 = date2jd(1950);
    jd1985  = date2jd(1985);
    daydiff = jd1950 - jd1985;
    tsec85 = (double(tday50) - daydiff)*86400;
    
end

