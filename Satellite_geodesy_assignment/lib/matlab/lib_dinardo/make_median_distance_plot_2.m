
function make_median_distance_plot_2(dist,RDS,SAR,MOD,type,unit,Mission_ID,Regional_Model_ID,RDS_LABEL)

numRDSD2C=sum(not(isnan(RDS)));
numSARD2C=sum(not(isnan(SAR)));
numMODD2C=sum(not(isnan(MOD)));

DIST=[150:200:9999];
best_fit=zeros(size(DIST));
for j=1:length(DIST)
    best_fit(j)=nanmedian(SAR(dist>DIST(j)-200 & dist<DIST(j)+200));       
end


bbest_fit=zeros(size(DIST));
for j=1:length(DIST)
    bbest_fit(j)=nanmedian(RDS(dist>DIST(j)-200 & dist<DIST(j)+200));       
end

bbbest_fit=zeros(size(DIST));

for j=1:length(DIST)
    
    
    bbbest_fit(j)=nanmedian(MOD(dist>DIST(j)-200 & dist<DIST(j)+200));
    
    
end



figure('Position',get(0,'Screensize'));hold on; 
%plot(dist,RDS,'b+');%plot PLRM
%hold on; plot(dist,SAR,'ko');% plot SAR
hold on;grid on
plot(DIST,bbest_fit,'g','LineWidth',3,'MarkerSize',8);
plot(DIST,best_fit,'r','LineWidth',3);
plot(DIST,bbbest_fit,'y','LineWidth',3,'MarkerSize',8);
set(gca,'FontSize',20);xlabel('DISTANCE [m]','FontSize',20);
ylabel([ type ' [' unit ']'],'FontSize',20);axis tight; 
%hl=legend ([ Mission_ID ' PLRM' ] ,[ Mission_ID ' SAR'],[ Mission_ID ' PLRM Median Curve'], [Mission_ID ' SAR Median Curve'],[ Regional_Model_ID ' Model Median Curve']);
hl=legend ([ Mission_ID ' ' RDS_LABEL ' Median Curve'], [Mission_ID ' SAR GPOD Median Curve'],[ Regional_Model_ID ' Model Median Curve']);
set(gca,'FontSize',20);set(hl,'Interpreter','none','FontSize',16);
text(0.2,0.90,['NP: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% text(0.2,0.90,['NP SAR: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% text(0.2,0.85,['NP PLRM: ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% text(0.2,0.80,['NP MODEL: ' num2str(numMODD2C,'% 10.3d')  ],'units','normalized','FontSize',22);

