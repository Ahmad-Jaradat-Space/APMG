
function make_percentile_plot(dist,RDS,SAR,type,unit,Mission_ID,RDS_LABEL)

numRDSD2C=sum(not(isnan(RDS)));
numSARD2C=sum(not(isnan(SAR)));

DIST=[150:200:9999];
prc_fit_rds=zeros(size(DIST));

for j=1:length(DIST)
    
    dummy=(RDS(dist>DIST(j)-200 & dist<DIST(j)+200));     
    npoint=sum(abs(dummy)<1);
    prc_fit_rds(j)=npoint/length(dummy)*100;
end



prc_fit_sar=zeros(size(DIST));

for j=1:length(DIST)
    
    dummy=(SAR(dist>DIST(j)-200 & dist<DIST(j)+200));     
    npoint=sum(abs(dummy)<1);
    prc_fit_sar(j)=npoint/length(dummy)*100;
end


figure('Position',get(0,'Screensize'));hold on; 

hold on; h=plot(DIST,prc_fit_sar,'r-','LineWidth',3,'MarkerSize',8);% plot SAR
hold on;grid on


plot(DIST,prc_fit_rds,'g-','LineWidth',3,'MarkerSize',8);

set(gca,'FontSize',20);xlabel('DISTANCE [m]','FontSize',16);
ylabel([ type ' [' unit ']'],'FontSize',20);
hl=legend ([ Mission_ID ' SAR GPOD'],[ Mission_ID ' ' RDS_LABEL ' ' ]  ,[ Mission_ID ' ' RDS_LABEL ' Median Curve'] ,[ Mission_ID ' SAR GPOD Median Curve']);set(gca,'FontSize',16);

text(0.4,0.90,['NP : ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);set(hl,'Interpreter','none');
% text(0.4,0.90,['NP SAR: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);set(hl,'Interpreter','none');
% text(0.4,0.85,['NP RDS: ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',22);

