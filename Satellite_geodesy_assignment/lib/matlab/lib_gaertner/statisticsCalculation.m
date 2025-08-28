function [stat,posSum,posHuncSum,dataStruct1,dataStruct2] = statisticsCalculation(dataStruct1,dataStruct2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


idrads = zeros(length(dataStruct1.LAT),1);
idbmle3 = zeros(length(dataStruct1.LAT),1);
for j=1:length(dataStruct1.LAT)
    [dummy, idx]=min(abs(dataStruct1.LAT(j)-dataStruct2.LAT));
    if dummy < 0.0005
        idbmle3(j)=idx;
        idrads(j)=j;
    end
end
idbmle3(idbmle3==0)=[] ;
idrads(idrads==0)=[] ;
dataStruct1.HUNC=dataStruct1.HUNC(idrads);
dataStruct2.HUNC=dataStruct2.HUNC(idbmle3);
pos = ~isnan(dataStruct1.HUNC) & ~isnan(dataStruct2.HUNC);
diff_hunc = dataStruct1.HUNC - dataStruct2.HUNC;
pos_hunc = pos;
dataStruct1.SWH=dataStruct1.SWH(idrads);
dataStruct2.SWH=dataStruct2.SWH(idbmle3);
diff_swh = dataStruct1.SWH - dataStruct2.SWH;
pos_swh = pos;
dataStruct1.U10=dataStruct1.U10(idrads);
dataStruct2.U10=dataStruct2.U10(idbmle3);
diff_u10 = dataStruct1.U10 - dataStruct2.U10;
pos_u10 = pos;
dataStruct1.sig0=dataStruct1.sig0(idrads);
dataStruct2.sig0=dataStruct2.sig0(idbmle3);
diff_sig0 = dataStruct1.sig0 - dataStruct2.sig0;
pos_sig0 = pos;
dataStruct1.LAT=dataStruct1.LAT(idrads);
dataStruct2.LAT=dataStruct2.LAT(idbmle3);

stat.bias.h_unc = nanmean(diff_hunc(pos_hunc));
stat.bias.swh = nanmean(diff_swh(pos_swh));
stat.bias.u10 = nanmean(diff_u10(pos_u10));
stat.bias.sig0 = nanmean(diff_sig0(pos_sig0));

stat.stdd.h_unc = nanstd(diff_hunc(pos_hunc));
stat.stdd.swh = nanstd(diff_swh(pos_swh));
stat.stdd.u10 = nanstd(diff_u10(pos_u10));
stat.stdd.sig0 = nanstd(diff_sig0(pos_sig0));

stat.corr.h_unc = corr(dataStruct1.HUNC(pos_hunc),dataStruct2.HUNC(pos_hunc));
% error here
stat.corr.swh = corr(dataStruct1.SWH(pos_swh),dataStruct2.SWH(pos_swh));
stat.corr.u10 = corr(dataStruct1.U10(pos_u10),dataStruct2.U10(pos_u10));
stat.corr.sig0 = corr(dataStruct1.sig0(pos_sig0),dataStruct2.sig0(pos_sig0));

posSum=sum(pos);
posHuncSum=sum(pos_hunc);

end

