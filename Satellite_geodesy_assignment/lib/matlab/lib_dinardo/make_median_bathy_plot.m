
function make_median_bathy_plot(bat,RDS,SAR,type,unit,Mission_ID,RDS_LABEL)

numRDSD2C=sum(not(isnan(RDS)));
numSARD2C=sum(not(isnan(SAR)));

BAT=[0:5:400];
best_fit=zeros(size(BAT));

for j=1:length( BAT)
    
    
    best_fit(j)=nanmedian(SAR(bat>BAT(j)-2.5 & bat<BAT(j)+2.5));
    
    
end

bbest_fit=zeros(size(BAT));

for j=1:length( BAT)
    
    
    bbest_fit(j)=nanmedian(RDS(bat>BAT(j)-2.5 & bat<BAT(j)+2.5));
    
    
end

figure('Position',get(0,'Screensize'));hold on; 
plot(bat,RDS,'b+');%plot PLRM
hold on; plot(bat,SAR,'ko');% plot SAR
hold on;grid on
plot(BAT,bbest_fit,'g-','LineWidth',3,'MarkerSize',8);
plot(BAT,best_fit,'r-','LineWidth',3,'MarkerSize',8);
set(gca,'FontSize',20);xlabel('BATHYMETRY [m]','FontSize',20);
ylabel([ type ' [' unit ']'],'FontSize',20);axis tight; 
hl=legend ([Mission_ID ' ' RDS_LABEL], [Mission_ID ' SAR GPOD' ],[Mission_ID ' ' RDS_LABEL ' Median Curve' ],[Mission_ID ' SAR GPOD Median Curve']);set(gca,'FontSize',20);
text(0.4,0.90,['NP SAR GPOD: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
text(0.4,0.85,['NP ' RDS_LABEL ': ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',22);set(hl,'Interpreter','none');

