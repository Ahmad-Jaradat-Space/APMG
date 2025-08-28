function [tsec] = mjd_to_sec1985(tmjd)
% 
% convert time from mjd to sec1985 
% 18.04.2014 L. Fenoglio new
%
    jdmjd = date2jd(1858, 11, 17);
    jd1985  = date2jd(1985);
    tsec = (tmjd - jd1985 + jdmjd)*86400.;
end

