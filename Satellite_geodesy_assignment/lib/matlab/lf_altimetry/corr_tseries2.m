function [rho,ave,stdval,minval,maxval]=corr_tseries2(x1,y1,x2,y2)
% subroutine to % compute correlation between 2 time series
%i x1,y1 values 1st timeseries 
%i x2,y2 values 2nd time series 
%o corr     correlation
%o rmsval   rms
%o minval   min  
%o maxval  max
%diff = x1 ~= x2;
diff=x1-x2;
rho=nancor(y1,y2);
diffalti=y1-y2;
ave=nanmean(diffalti);
stdval=nanstd(diffalti);
maxval=nanmax(diffalti);
minval=nanmin(diffalti);