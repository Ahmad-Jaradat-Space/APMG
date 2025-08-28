
function make_std_bathy_plot_2(bat,RDS,SAR,MOD,type,unit,Mission_ID,Regional_Model_ID,RDS_LABEL)

numRDSD2C=sum(not(isnan(RDS)));
numSARD2C=sum(not(isnan(SAR)));
numMODD2C=sum(not(isnan(MOD)));

BAT=[0:1:100];
best_fit=zeros(size(BAT));
for j=1:length(BAT)
    best_fit(j)=nanstd(SAR(bat>BAT(j)-0.5 & bat<BAT(j)+0.5));       
end


bbest_fit=zeros(size(BAT));
for j=1:length(BAT)
    bbest_fit(j)=nanstd(RDS(bat>BAT(j)-0.5 & bat<BAT(j)+0.5));       
end

bbbest_fit=zeros(size(BAT));

for j=1:length(BAT)
    
    
    bbbest_fit(j)=nanstd(MOD(bat>BAT(j)-0.5 & bat<BAT(j)+0.5));
    
    
end



figure('Position',get(0,'Screensize'));hold on; 
hold on;grid on
bar(BAT,bbest_fit,'histc');
bar(BAT,best_fit,'histc');
bar(BAT,bbbest_fit,'histc')
h = findobj(gca,'Type','patch');
set(gcf,'Renderer','opengl')
%set(h(1),'EdgeColor','k');hold on;
set(h(1),'FaceColor','g');set(gca,'FontSize',20)
set(h(2),'FaceColor','r');set(gca,'FontSize',20)
set(h(3),'FaceColor','b');set(gca,'FontSize',20)
set(gca,'FontSize',20);xlabel('BATHYMETRY [m]','FontSize',20);
ylabel([ type ' [' unit ']'],'FontSize',20);axis tight; 
hl=legend ([ Mission_ID ' PLRM'],[Mission_ID ' SAR GPOD'],[ Regional_Model_ID ' Ocean Model']);
set(hl,'Interpreter','none');
set(gca,'FontSize',20);axis tight
text(0.4,0.90,['NP SAR GPOD: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
text(0.4,0.85,['NP ' RDS_LABEL ': ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
text(0.4,0.80,['NP MODEL: ' num2str(numMODD2C,'% 10.3d')  ],'units','normalized','FontSize',22);

