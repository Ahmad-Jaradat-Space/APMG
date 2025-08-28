
function make_std_distance_plot_2(dist,RDS,SAR,MOD,type,unit,Mission_ID,Regional_Model_ID,RDS_LABEL)

numRDSD2C=sum(not(isnan(RDS)));
numSARD2C=sum(not(isnan(SAR)));
numMODD2C=sum(not(isnan(MOD)));

DIST=[150:200:9999];
best_fit=zeros(size(DIST));
for j=1:length(DIST)
    best_fit(j)=nanstd(SAR(dist>DIST(j)-200 & dist<DIST(j)+200));       
end


bbest_fit=zeros(size(DIST));
for j=1:length(DIST)
    bbest_fit(j)=nanstd(RDS(dist>DIST(j)-200 & dist<DIST(j)+200));       
end

bbbest_fit=zeros(size(DIST));

for j=1:length(DIST)
    
    
    bbbest_fit(j)=nanstd(MOD(dist>DIST(j)-200 & dist<DIST(j)+200));
    
    
end



figure('Position',get(0,'Screensize'));hold on; 
hold on;grid on
% 
% plot(DIST,bbest_fit,'g');
% plot(DIST,best_fit,'r');
% plot(DIST,bbbest_fit,'b')

bar(DIST,bbest_fit,'histc');
bar(DIST,best_fit,'histc');
bar(DIST,bbbest_fit,'histc')
h = findobj(gca,'Type','patch');
set(gcf,'Renderer','opengl')
%set(h(1),'EdgeColor','k');hold on;
set(h(1),'FaceColor','g');set(gca,'FontSize',20)
set(h(2),'FaceColor','r');set(gca,'FontSize',20)
set(h(3),'FaceColor','b');set(gca,'FontSize',20)

set(gca,'FontSize',20);xlabel('DISTANCE [m]','FontSize',20);
ylabel([ type ' [' unit ']'],'FontSize',20);axis tight; 

if not(nansum(MOD)==0)
    
hl=legend ([ Mission_ID ' ' RDS_LABEL ], [ Mission_ID ' SAR GPOD'],[Regional_Model_ID ' Model' ]);

else
    
    hl=legend ([ Mission_ID ' ' RDS_LABEL ], [ Mission_ID ' SAR GPOD']);

end

set(hl,'Interpreter','none');

% if not(nansum(MOD)==0)
% 
% set(gca,'FontSize',20);axis tight
% text(0.4,0.90,['NP SAR: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% text(0.4,0.85,['NP RDS: ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% text(0.4,0.80,['NP MODEL: ' num2str(numMODD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% 
% else
%     
%  set(gca,'FontSize',20);axis tight
% text(0.4,0.90,['NP SAR: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
% text(0.4,0.85,['NP RDS: ' num2str(numRDSD2C,'% 10.3d')  ],'units','normalized','FontSize',22);
%  
% end


text(0.8,0.70,['NP: ' num2str(numSARD2C,'% 10.3d')  ],'units','normalized','FontSize',22);


