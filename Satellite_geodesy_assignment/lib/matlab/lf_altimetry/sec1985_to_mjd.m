function [tmjd] = sec1985_to_mjd(tsec)
% Transfor sec1985 to mjd 
% 18.04.2014 L. Fenoglio new
%  
% convert time from sec1985 to mjd 
%
    jdmjd = date2jd(1858, 11, 17);
    jd1985  = date2jd(1985);
    tmjd = (tsec/86400. + jd1985 - jdmjd);
end

