function [xmov,ymov,xia,yia,trendmov,trendia]=ts2mov(x1,y1,lag)
% subroutine to % compute correlation between 2 time series
% given a time series, do the following steps:
% - demean,
% - deseasonalize,
% running mean
%
% demean
flag_filtre=0;
y1=y1-nanmean(y1);
lag=12;
% transform from column to raw vectors date and ts !
date=x1';
ts=y1';
% eliminate the annual and semiannual signal from altimetry msl
if (flag_filtre == 0)
 [ts_out] = ts_extr_an_seman(date,ts,2);
else
 [ts_out] = ts_filtre_hf(date,ts,1/12,fref);
end
mslts=ts_out;
mslts2=mslts-mean(mslts);
% compute trend of the ia tseries
%
[ts_detrend,trendia] = ts_extr_trend(x1',mslts2);
%trendia=0.d0;
% compute running average over the residuals
xia = x1;
yia=mslts2;
[nbl,nbc] = size(yia);
[ymov]=moy_mouv(yia,lag);
xmov = xia(lag/2:(end-lag/2));
%trendmov=0.d0;
% compute trend of the running average
[ts_detrend,trendmov] = ts_extr_trend(xmov',ymov);
% trendmov