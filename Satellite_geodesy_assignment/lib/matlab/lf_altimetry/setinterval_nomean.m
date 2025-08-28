function [x1,y1]=setinterval_nomean(x0,y0,t1,t2)
% subroutine to extract a subinterval  
%i x0 time series
%i y0 corresponding values
%i t1 start time
%i t2 end time
%o x1 time in output 
%o y1 values in output
datey=x0;
posalti=find(datey> t1 & datey<t2);
x1=x0(posalti);
y1=y0(posalti);
y1 = y1 ;