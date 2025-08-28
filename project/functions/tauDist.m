function tau = tauDist(alpha,r,n,p)

% Evaluates the TAU distribution
%
% Input:
% alpha             scalar, double, significance level
% r                 scalar, int, redundancy
% n                 scalar, int, number of observations
% p                 scalar, int, number of tests performed (default = 1)
%
% Output:
% tau               scalar, double, critical value of the TAU distribution
%
%%       References:
%       [1] Koch KR,1999, Parameter estimation and hypothesis testing in linear models,
%       Springer Verlag, Berlin
%       [2] Pope AJ,1976, The  statistics  of residuals and the detection of outliers.
%       NOAA Technical Report NOS65 NGS1, US Department of Commerce, National Geodetic
%       Survey,Rockville, Maryland
%
% Formulae from Koch (1999) 4.117

if ~exist('p','var')                   % p given?
    p=1;
end
% Fisher distribution percentile
fdist = finv((1-alpha/n),p,(r-p));
% TAU distribution evaluation
tau = sqrt((r*fdist)/(r-p+p*fdist));
