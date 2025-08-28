function plotting(config,dataStruct1,dataStruct2,stat,posSum,posHuncSum)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


figure('Position',get(0,'Screensize'));
subplot(4,1,1)
plot(dataStruct2.LAT,dataStruct2.HUNC,'Marker','.','Markeredgecolor','b','Linestyle','none');
hold on
plot(dataStruct1.LAT,dataStruct1.HUNC,'Marker','.','Markeredgecolor','r','Linestyle','none');
axis tight
xlabel('Latitude [deg]','Fontsize',16)
ylabel('SSH_{unc} [m]','Fontsize',16)
legend(config.nameFile1,config.nameFile2)
text(0.8,0.9,sprintf('Bias = %.3fm',stat.bias.h_unc),'Units','normalized','Fontsize',16)
text(0.8,0.75,sprintf('Stdd = %.3fm',stat.stdd.h_unc),'Units','normalized','Fontsize',16)
text(0.8,0.6,sprintf('Corr = %.3f',stat.corr.h_unc),'Units','normalized','Fontsize',16)
text(0.8,0.45,sprintf('N_{start}/N_{end} = %d/%d',posSum,posHuncSum),'Units','normalized','Fontsize',16)
set(gca,'Fontsize',16)
hold on


subplot(4,1,2)
plot(dataStruct2.LAT,dataStruct2.SWH,'Marker','.','Markeredgecolor','b','Linestyle','none');
hold on
plot(dataStruct1.LAT,dataStruct1.SWH,'Marker','.','Markeredgecolor','r','Linestyle','none');
axis tight
xlabel('Latitude [deg]','Fontsize',16)
ylabel('SWH [m]','Fontsize',16)
legend(config.nameFile1,config.nameFile2)
text(0.8,0.9,sprintf('Bias = %.3fm',stat.bias.swh),'Units','normalized','Fontsize',16)
text(0.8,0.75,sprintf('Stdd = %.3fm',stat.stdd.swh),'Units','normalized','Fontsize',16)
text(0.8,0.6,sprintf('Corr = %.3f',stat.corr.swh),'Units','normalized','Fontsize',16)
text(0.8,0.45,sprintf('N_{start}/N_{end} = %d/%d',posSum,posHuncSum),'Units','normalized','Fontsize',16)
set(gca,'Fontsize',16)
hold on


subplot(4,1,3)
plot(dataStruct2.LAT,dataStruct2.U10,'Marker','.','Markeredgecolor','b','Linestyle','none');
hold on
plot(dataStruct1.LAT,dataStruct1.U10,'Marker','.','Markeredgecolor','r','Linestyle','none');
axis tight
xlabel('Latitude [deg]','Fontsize',16)
ylabel('U10 [m/s]','Fontsize',16)
legend(config.nameFile1,config.nameFile2)
text(0.8,0.9,sprintf('Bias = %.3fm/s',stat.bias.u10),'Units','normalized','Fontsize',16)
text(0.8,0.75,sprintf('Stdd = %.3fm/s',stat.stdd.u10),'Units','normalized','Fontsize',16)
text(0.8,0.6,sprintf('Corr = %.3f',stat.corr.u10),'Units','normalized','Fontsize',16)
text(0.8,0.45,sprintf('N_{start}/N_{end} = %d/%d',posSum,posHuncSum),'Units','normalized','Fontsize',16)
set(gca,'Fontsize',16)
hold on


subplot(4,1,4)
plot(dataStruct2.LAT,dataStruct2.sig0,'Marker','.','Markeredgecolor','b','Linestyle','none');
hold on
plot(dataStruct1.LAT,dataStruct1.sig0,'Marker','.','Markeredgecolor','r','Linestyle','none');
axis tight
xlabel('Latitude [deg]','Fontsize',16)
ylabel('sig0 [dB]','Fontsize',16)
legend(config.nameFile1,config.nameFile2)
text(0.8,0.9,sprintf('Bias = %.3fdB',stat.bias.sig0),'Units','normalized','Fontsize',16)
text(0.8,0.75,sprintf('Stdd = %.3fdB',stat.stdd.sig0),'Units','normalized','Fontsize',16)
text(0.8,0.6,sprintf('Corr = %.3f',stat.corr.sig0),'Units','normalized','Fontsize',16)
text(0.8,0.45,sprintf('N_{start}/N_{end} = %d/%d',posSum,posHuncSum),'Units','normalized','Fontsize',16)
set(gca,'Fontsize',16)
saveas(gcf,['res/AlongTrack/test_3Par_',config.nameFile1,'_',config.nameFile2,'.png'])
hold on




figure;
plot(dataStruct1.LON,dataStruct1.LAT,'Marker','.','Markeredgecolor','r','Linestyle','none');
axis tight
xlabel('Lon [deg]','Fontsize',16)
ylabel('Lat [deg]','Fontsize',16)
legend('S3A test track')
set(gca,'Fontsize',16)
saveas(gcf,['res/AlongTrack/test_Track_',config.nameFile1,'_',config.nameFile2,'.png'])


end

