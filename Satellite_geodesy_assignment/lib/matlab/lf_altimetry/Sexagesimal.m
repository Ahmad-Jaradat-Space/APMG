%SEXAGESIMAL Converts decimal to sexagesimal number
%
% INPUT
%   dez - decimal number (deg)
%
% OUTPUT
%   sgz - sexagesimal number (deg.arcminArcsec)
%
% EXAMPLE
%   sgz = Sexagesimal(8.50639);
%   sgz = 8.3023;
%
% AUTHOR: F. Reichel
% DATE: 01.04.2015
function sgz = Sexagesimal(dez)
sgz = fix(dez);
frac = dez-sgz;
am = fix(frac*60);
as = (frac*60-am)*60;
sgz = sgz+am/100+as/10000;