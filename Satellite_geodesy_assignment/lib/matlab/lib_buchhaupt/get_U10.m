function [U10,std_U10] = get_U10(sigma0,KA,std_sigma0)
% This is a 1-D wind speed model relating backscatter (sigma0) to wind speed at
% 10 meters altitude. It was developed at ECMWF by Saleh Abdalla to support the
% Envisat mission.
%
% To support also the SARAL processing, the optional parameter <ka_band> was added,
% which will select parameters appropriate for Ka-band instead of Ku-band.
%
% References:
% Abdalla, S., Ku-band radar altimeter surface wind speed algorithm,
%  in Proc. of the 2007 Envisat Symposium, Montreux, Switzerland, 23-27 April 2007,
%  Eur. Space Agency Spec. Publ., ESA SP-636, 2007.
% Lillibridge, J. L., R. Scharroo, S. Abdalla, and D. C. Vandemark, One- and
%  two-dimensional wind speed models for Ka-band altimetry, J. Atmos. Oceanic
%  Technol., 31(3), 630-638, doi:10.1175/JTECH-D-13- 00167.1, 2014.
%
% Input:
%  sigma0 : backscatter coefficient (dB)
%  KA: backscatter is from Ka-band instead of Ku-band
%
% Output:
%  U10 : wind speed at 10 meters altitude (m/s)

if nargin < 2
    KA = false;
end

if KA % Ka-band values (c and sigmab modified for continuity)
    a = 34.2;
    b = 2.48;
    c = 711.6;
    d = 0.42;
    sigmab = 11.409;
else % Ku-band values
    a = 46.5;
    b = 3.6;
    c = 1690;
    d = 0.5;
    sigmab = 10.917;
end

U10 = nan(size(sigma0));
std_U10 = U10;
um = U10;

% Step 1: First approximation
for i=1:length(sigma0)
    if isnan(sigma0(i))
        % If SIGMA0 is NaN, return NaN
        U10(i) = sigma0(i);
        if nargout == 2
            std_U10(i) = nan;
        end
    elseif (sigma0(i) <= sigmab)
        % Linear portion of the model
        um(i) = a - b * sigma0(i);
        if nargout == 2
            std_U10(i) = b*std_sigma0(i);
        end
    else
        % Exponential portion of the model
        um(i) = c * exp(-d * sigma0(i));
        if nargout == 2
            std_U10(i) = c*d*std_sigma0(i)*exp(-d*sigma0(i));
        end
    end
end

% Step 2: Fine tuning

% ux = um^0.096;
ux = um.^0.096;
U10 = um + 1.4 * ux .* exp(-0.32 .* um .* ux);

end

