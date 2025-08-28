function mjd = decyear2mjd(decyear)
% decyear2mjd - Convert decimal year to Modified Julian Date
%
% SYNTAX:
%   mjd = decyear2mjd(decyear)
%
% INPUT:
%   decyear - Decimal year (e.g., 2005.5 for middle of 2005)
%
% OUTPUT:
%   mjd     - Modified Julian Date
%
% NOTES:
%   - MJD = JD - 2400000.5
%   - Reference epoch: January 1, 1858, 00:00 UTC
%   - Accounts for leap years
%
% EXAMPLE:
%   mjd = decyear2mjd(2005.0)  % January 1, 2005
%   mjd = decyear2mjd(2005.5)  % Middle of 2005
%
% Author: GRACE Analysis Project
% Date: 2025

% Input validation
if nargin < 1
    error('Decimal year input required');
end

% Handle vector inputs
if length(decyear) > 1
    mjd = zeros(size(decyear));
    for i = 1:length(decyear)
        mjd(i) = decyear2mjd(decyear(i));
    end
    return;
end

% Extract integer year and fractional part
year = floor(decyear);
frac_year = decyear - year;

% Check if it's a leap year
is_leap = (mod(year, 4) == 0 & mod(year, 100) ~= 0) | (mod(year, 400) == 0);

% Days in year
if is_leap
    days_in_year = 366;
else
    days_in_year = 365;
end

% Convert to MATLAB datenum (days since January 0, 0000)
jan1_datenum = datenum(year, 1, 1);

% Add fractional year contribution
total_datenum = jan1_datenum + frac_year * days_in_year;

% Convert to MJD
% MJD reference: November 17, 1858 (MJD = 0)
% MATLAB datenum reference: January 1, 0000 (datenum = 1)
mjd_reference_datenum = datenum(1858, 11, 17);
mjd = total_datenum - mjd_reference_datenum;

end