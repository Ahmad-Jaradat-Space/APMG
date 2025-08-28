%DEZIMAL Converts sexagesimal to decimal number
%
% INPUT
%   sgz - sexagesimal number (deg.arcminArcsec)
%
% OUTPUT
%   dez - decimal number (deg)
%
% EXAMPLE
%   dez = Dezimal(8.3023);
%   dez = 8.50639;
%
% AUTHOR: F. Reichel
% DATE: 01.04.2015
function dez = Dezimal(sgz)
dez = fix(sgz);
frac = sgz - dez;
am = fix(frac * 100);
as = (frac * 100 - am) * 100;
dez = dez + am/60 + as/3600;