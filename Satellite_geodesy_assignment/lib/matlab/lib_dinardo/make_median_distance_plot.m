
function make_median_distance_plot(dist,RDS,SAR,type,unit,Mission_ID,RDS_LABEL)

numRDSD2C=sum(not(isnan(RDS)));
numSARD2C=sum(not(isnan(SAR)));

DIST=[150:200:9999];
best_fit=zeros(size(DIST));
for j=1:length(DIST)
    best_fit(j)=nanmedian(SAR(dist>DIST(j)-200 & dist<DIST(j)+200));       
end


bbest_fit=zeros(size(DIST));
for j=1:length(DIST)
    bbest_fit(j)=nanmedian(RDS(dist>DIST(j)-200 & dist<DIST(j)+200));       
end

figure('Position',get(0,'Screensize'));hold on; 

h1 = subplot(2,1,1);

plot(dist,RDS,'bo');%plot PLRM
hold on;grid on
%plot(DIST,bbest_fit,'gs-','LineWidth',3,'MarkerSize',8);
plot(DIST,bbest_fit,'g-','LineWidth',3,'MarkerSize',8);
plot(DIST,best_fit,'r','LineWidth',3,'MarkerSize',8);
set(gca,'FontSize',20);xlabel('DISTANCE [m]','FontSize',16);
ylabel([ type ' [' unit ']'],'FontSize',20);axis tight; 
hl=legend (h1, [ Mission_ID ' ' RDS_LABEL ' Points' ] ,[ Mission_ID ' ' RDS_LABEL ' Median Curve'] ,[ Mission_ID ' SAR GPOD Median Curve']);set(gca,'FontSize',16);
%text(0.45,0.90,['NP SAR: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',18);set(hl,'Interpreter','none');
%text(0.45,0.8,['NP PLRM: ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',18);
text(0.45,0.90,['NP: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',18);set(hl,'Interpreter','none');

h2 = subplot(2,1,2);

plot(dist,SAR,'ko');% plot SAR;%plot PLRM
%hold on; plot(dist,SAR,'ko');% plot SAR
hold on;grid on
%plot(DIST,bbest_fit,'gs-','LineWidth',3,'MarkerSize',8);
%plot(DIST,bbest_fit,'g-','LineWidth',3,'MarkerSize',8);
plot(DIST,best_fit,'r','LineWidth',3,'MarkerSize',8);
set(gca,'FontSize',20);xlabel('DISTANCE [m]','FontSize',16);


ylabel([ type ' [' unit ']'],'FontSize',20);axis tight; 
hl=legend (h2, [ Mission_ID ' SAR GPOD Points'],[ Mission_ID ' ' RDS_LABEL  ' Median Curve'] ,[ Mission_ID ' SAR GPOD Median Curve']);set(gca,'FontSize',16);
text(0.45,0.90,['NP: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',18);set(hl,'Interpreter','none');
%text(0.45,0.90,['NP SAR: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',18);set(hl,'Interpreter','none');
%text(0.45,0.8,['NP PLRM: ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',18);

set(gcf, 'currentaxes', h1)