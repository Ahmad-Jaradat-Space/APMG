
function [tt hh] = fillNaN(t,h,dt)

%% PURPOSE:
% fill data gaps of regular-spaced time series with NaN

%% INPUT:
% t: time info of exsisting time series [days]
% h: data info of exsisting time series
% dt: time interval of time series

%% OUTPUT:
% tt: time steps of new time series
% hh: data of new time series

%%

% extend of input time series:
tt=(min(t):dt:max(t))';  

% number of regular points:
no=size(tt,1);

% initial settings: all time series elements =NaN
hh=ones(no,1);
hh(:)=NaN;

% fill new time series with exsisting data
if (no>size(t,1))
    for i=1:no
        index=find(t==tt(i));     
        if (index)
           hh(i)=h(index);
        end
    end
end