%UTM_GK2BL Transforms UTM or GaussKrueger coordinates to ellipsoidal latitude,longitude (BL)
% INPUT
%   X     - UTM North or GK Hoch coordinate (m)
%   Y     - UTM East or GK Rechts coordinate (m)
%   plane - 'UTM' or 'GK'
%
% OUTPUT
%   BL - 2x1 Vector (B;L) (sexagesimal, sexagesimal)
%
% EXAMPLE
%   BL = UTM_GK2BL(5523330.7658,32476901.1292,'UTM')
%
% AUTHOR: F. Reichel
% DATE: 01.04.2015

function BL = UTM_GK2BL(X,Y,plane)

if nargin < 3
    exit;
end

if strcmp(plane,'UTM')
    ellipsoid='GRS80';
    m=0.9996;
    centMeridID=fix(Y/1E6);                 % central meridian zone number 
    sgzL0=Sexagesimal(centMeridID*6-183);   % central meridian
elseif strcmp(plane,'GK')
    ellipsoid='Bessel';
    m=1;
    centMeridID=fix(Y/1E6);                 % central meridian zone number
    sgzL0=Sexagesimal(centMeridID*3);       % central meridian
else
    disp('Plane unknown! Use UTM or GK.');
    exit
end

X = X/m;
Y = (Y - (centMeridID + 0.5) * 1E6)/m;

BL = XY2BL(X,Y,sgzL0,ellipsoid);